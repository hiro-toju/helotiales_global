#!/usr/bin/env python3
"""Reproducible ITS/rRNA phylogeny and phylogenetic-signal workflow."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback
from collections import Counter, defaultdict
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Iterable, Iterator, Optional, TextIO

import pandas as pd


DEFAULT_TRAIT = "host_z_standardized_dprime_sampling_unit_label_shuffle"
PLATFORM_ORDER = {
    "Sanger": 0, "PacBio": 1, "Illumina": 2, "DNBSEQ": 3,
    "Ion Torrent": 4, "454 Roche": 5, "Other": 6, "Unknown": 7,
}
MISSING_TAXA = {"", "unknown", "unidentified", "uncultured", "fungi sp", "fungi_sp"}
RUN_STARTED_AT = datetime.now(timezone.utc).astimezone()
RUN_STARTED_MONOTONIC = time.perf_counter()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, default=Path("data"))
    p.add_argument("--summary", type=Path)
    p.add_argument("--two-dp-table", type=Path)
    p.add_argument("--input-fasta", type=Path)
    p.add_argument("--unite", type=Path)
    p.add_argument(
        "--output-root", type=Path, default=Path("outputs"),
        help=(
            "Root directory for automatically named outputs. When --specificity-run-name "
            "is set, the actual output directory is <output-root>/<run-name>."
        ),
    )
    p.add_argument(
        "--specificity-run-name",
        required=True,
        help=(
            "Run folder name produced by a specificity workflow, for example "
            "sh150_family30_10000 or sh150_family30_10000_lat20_minblock100. "
            "Host-specificity inputs and FASTA are taken from that run."
        ),
    )
    p.add_argument(
        "--specificity-output-root", type=Path,
        default=Path("../4_Specificity_Lat/outputs"),
        help=(
            "Root directory containing specificity output run folders. "
            "Default: ../4_Specificity_Lat/outputs. Use "
            "../3_Specificity_SpatialBlocks/outputs to analyze SpatialBlocks results."
        ),
    )
    p.add_argument(
        "--specificity-fasta-subset",
        choices=("its1_and_its2", "all", "its1_only", "its2_only", "selected"),
        default="its1_and_its2",
        help=(
            "FASTA subset to use from the specificity run. "
            "The default its1_and_its2 uses the ITS1+ITS2 FASTA copied from 1_GlobalFungi. "
            "Use selected to reproduce the run-specific selected_SH_sequences.fasta, "
            "which may include fallback sequences from other ITS subsets."
        ),
    )
    p.add_argument("--trait-column", default=DEFAULT_TRAIT)
    p.add_argument("--target-taxon", default="Helotiales")
    p.add_argument("--minimum-length", type=int, default=400)
    p.add_argument(
        "--outgroup-minimum-length", type=int, default=460,
        help=(
            "Minimum ungapped A/C/G/T length for UNITE outgroup candidates "
            "before BLAST and alignment (default: 460)"
        ),
    )
    p.add_argument(
        "--exclude-sequence-name-keyword", action="append", default=[], metavar="TEXT",
        help="Exclude target FASTA records whose header contains TEXT; repeatable",
    )
    p.add_argument("--disable-default-sequence-name-exclusions", action="store_true")
    p.add_argument("--sequence-name-keyword-case-sensitive", action="store_true")
    p.add_argument("--maximum-ambiguous-fraction", type=float, default=0.02)
    p.add_argument("--minimum-outgroup-query-coverage", type=float, default=50.0)
    p.add_argument("--outgroup-minimum-median-identity", type=float, default=75.0)
    p.add_argument("--outgroup-maximum-median-identity", type=float, default=85.0)
    p.add_argument("--outgroup-target-median-identity", type=float, default=80.0)
    p.add_argument("--outgroup-maximum-identity-to-any-target", type=float, default=90.0)
    p.add_argument("--outgroup-minimum-median-query-coverage", type=float, default=75.0)
    p.add_argument("--outgroup-minimum-query-hit-fraction", type=float, default=0.80)
    p.add_argument("--outgroup-required-shared-rank", choices=("class", "phylum"), default="class")
    p.add_argument("--minimum-outgroup-genus-sh-count", type=int, default=3)
    p.add_argument("--misidentification-minimum-pairwise-overlap", type=int, default=100)
    p.add_argument("--misidentification-minimum-comparisons", type=int, default=3)
    p.add_argument("--misidentification-mad-multiplier", type=float, default=3.0)
    p.add_argument("--target-vs-outgroup-mean-identity-margin", type=float, default=0.0)
    p.add_argument("--mafft-algorithm", choices=("linsi", "einsi", "ginsi", "auto"), default="linsi")
    p.add_argument(
        "--terminal-trim-minimum-occupancy", type=float, default=1.0,
        help=(
            "Minimum fraction of sequences with A/C/G/T required at the left and "
            "right terminal-trimming boundaries after MAFFT (default: 1.0)"
        ),
    )
    p.add_argument("--trimmer", choices=("trimal",), default="trimal")
    p.add_argument("--trimal-gap-threshold", type=float, default=0.50)
    p.add_argument(
        "--trimal-minimum-retained-percent", type=float, default=0.0,
        help=(
            "Minimum percentage of original alignment columns retained by trimAl "
            "(-cons); 0 disables this override (default)"
        ),
    )
    p.add_argument("--bootstrap", type=int, default=500)
    p.add_argument("--signal-randomizations", type=int, default=10000)
    p.add_argument("--nj-distance-model", default="K80")
    p.add_argument("--seed", type=int, default=20260705)
    p.add_argument("--jobs", type=int, default=0)
    p.add_argument("--iqtree-thread-mode", choices=("auto", "fixed"), default="auto")
    p.add_argument("--figure-font-scale", type=float, default=1.2)
    p.add_argument(
        "--iqtree-executable", type=Path,
        help="IQ-TREE executable; overrides project-local and PATH installations",
    )
    p.add_argument("--r-script", type=Path, default=Path(__file__).with_name("phylogeny_analysis.R"))
    args = p.parse_args()
    default_keywords = [] if args.disable_default_sequence_name_exclusions else ["Ascomycota_sp"]
    args.exclude_sequence_name_keyword = list(dict.fromkeys(
        [*default_keywords, *args.exclude_sequence_name_keyword]
    ))
    if args.minimum_length < 1 or args.outgroup_minimum_length < 1 or args.signal_randomizations < 1:
        p.error("length and signal-randomization values must be positive")
    if args.bootstrap < 100:
        p.error("--bootstrap must be at least 100 for IQ-TREE standard bootstrap")
    if not 0 <= args.maximum_ambiguous_fraction < 1:
        p.error("--maximum-ambiguous-fraction must be in [0, 1)")
    if not 0 <= args.minimum_outgroup_query_coverage <= 100:
        p.error("--minimum-outgroup-query-coverage must be in [0, 100]")
    if not (0 <= args.outgroup_minimum_median_identity
            < args.outgroup_target_median_identity
            < args.outgroup_maximum_median_identity <= 100):
        p.error("outgroup median-identity settings must satisfy 0 <= minimum < target < maximum <= 100")
    if not 0 <= args.outgroup_maximum_identity_to_any_target <= 100:
        p.error("--outgroup-maximum-identity-to-any-target must be in [0, 100]")
    if not 0 <= args.outgroup_minimum_median_query_coverage <= 100:
        p.error("--outgroup-minimum-median-query-coverage must be in [0, 100]")
    if not 0 < args.outgroup_minimum_query_hit_fraction <= 1:
        p.error("--outgroup-minimum-query-hit-fraction must be in (0, 1]")
    if args.minimum_outgroup_genus_sh_count < 1:
        p.error("--minimum-outgroup-genus-sh-count must be at least 1")
    if args.misidentification_minimum_pairwise_overlap < 1:
        p.error("--misidentification-minimum-pairwise-overlap must be at least 1")
    if args.misidentification_minimum_comparisons < 1:
        p.error("--misidentification-minimum-comparisons must be at least 1")
    if (not math.isfinite(args.misidentification_mad_multiplier)
            or args.misidentification_mad_multiplier <= 0):
        p.error("--misidentification-mad-multiplier must be finite and greater than 0")
    if (not math.isfinite(args.target_vs_outgroup_mean_identity_margin)
            or args.target_vs_outgroup_mean_identity_margin < 0):
        p.error("--target-vs-outgroup-mean-identity-margin must be finite and non-negative")
    if not 0 < args.terminal_trim_minimum_occupancy <= 1:
        p.error("--terminal-trim-minimum-occupancy must be in (0, 1]")
    if not 0 <= args.trimal_gap_threshold <= 1:
        p.error("--trimal-gap-threshold must be in [0, 1]")
    if not 0 <= args.trimal_minimum_retained_percent <= 100:
        p.error("--trimal-minimum-retained-percent must be in [0, 100]")
    if args.jobs < 0:
        p.error("--jobs must be 0 (automatic) or a positive integer")
    if not math.isfinite(args.figure_font_scale) or args.figure_font_scale <= 0:
        p.error("--figure-font-scale must be finite and greater than 0")
    return args


def one_match(directory: Path, pattern: str, label: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly one {label} matching {directory / pattern}; found {len(matches)}")
    return matches[0]


def resolve_existing_path(path: Path, project_dir: Path) -> Path:
    candidates = [path]
    if not path.is_absolute():
        candidates.extend([Path.cwd() / path, project_dir / path, project_dir.parent / path])
    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved.exists():
            return resolved
    return path.expanduser().resolve()


def find_run_directory(root: Path, run_name: str, label: str) -> Path:
    root = root.resolve()
    direct = root / run_name
    if direct.is_dir():
        return direct
    matches = sorted(path for path in root.glob(f"**/{run_name}") if path.is_dir())
    if len(matches) != 1:
        found = ", ".join(str(path) for path in matches[:8])
        raise FileNotFoundError(
            f"Expected exactly one {label} run directory named {run_name} under {root}; "
            f"found {len(matches)}. {found}"
        )
    return matches[0]


def resolve_specificity_inputs(args: argparse.Namespace, project_dir: Path) -> dict[str, Path]:
    run_name = args.specificity_run_name
    root = resolve_existing_path(args.specificity_output_root, project_dir)
    run = find_run_directory(root, run_name, "specificity")
    result_dir = one_match(
        run, "*root_occurrences_with_hosts",
        "specificity result directory",
    )
    if args.specificity_fasta_subset == "selected":
        input_fasta = one_match(
            result_dir, "*_by_host_*_selected_SH_sequences.fasta",
            "specificity selected SH FASTA",
        )
    else:
        input_fasta = one_match(
            run / "data",
            f"*root_host_sequences_{args.specificity_fasta_subset}.fasta",
            f"specificity copied {args.specificity_fasta_subset} FASTA",
        )
    return {
        "specificity_run_dir": run,
        "specificity_result_dir": result_dir,
        "summary": one_match(
            result_dir, "*_by_host_*_selected_fungi_integrated_summary.tsv",
            "specificity host integrated summary",
        ),
        "two_dp_table": one_match(
            result_dir, "*_by_host_*_sampling_unit_label_shuffle_2dp_z.tsv",
            "specificity host 2DP z table",
        ),
        "input_fasta": input_fasta,
    }


def clean_seq(value: str) -> str:
    return re.sub(r"[^ACGTRYSWKMBDHVN?-]", "", value.upper().replace("U", "T"))


def called_length(seq: str) -> int:
    return len(re.sub(r"[-?.]", "", seq))


def ambiguous_fraction(seq: str) -> float:
    called = re.sub(r"[-?.]", "", seq)
    if not called:
        return 1.0
    return sum(base not in "ACGT" for base in called) / len(called)


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header: Optional[str] = None
        chunks: list[str] = []
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, clean_seq("".join(chunks))
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line)
        if header is not None:
            yield header, clean_seq("".join(chunks))


def write_fasta(records: Iterable[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def safe_label(value: object) -> str:
    text = re.sub(r"\s+", "_", str(value).strip())
    return re.sub(r"[^A-Za-z0-9_.=+\-\[\]]", "_", text).strip("_") or "Unknown"


def sh_id(header: str) -> str:
    hit = re.search(r"SH\d+(?:\.\d+FU)?", header)
    return hit.group(0) if hit else header.split("|", 1)[0].split()[0]


def fasta_species(header: str) -> str:
    parts = header.split("|")
    return parts[1] if len(parts) > 1 else "Unknown"


def fasta_n(header: str) -> int:
    hit = re.search(r"(?:^|\|)n_unique_occurrences=(\d+)(?:\||$)", header)
    return int(hit.group(1)) if hit else 0


def command_version(executable: str) -> str:
    for flag in ("--version", "-version", "-h"):
        try:
            out = subprocess.run([executable, flag], capture_output=True, text=True, timeout=20)
            text = (out.stdout or out.stderr).strip().splitlines()
            if text:
                return text[0][:500]
        except Exception:
            pass
    return "version unavailable"


def portable_report_path(value: object, project_dir: Path) -> str:
    path = Path(value).expanduser()
    try:
        resolved = path.resolve()
    except OSError:
        return f"[EXTERNAL]/{path.name}"
    for base in (project_dir, project_dir.parent):
        try:
            return str(resolved.relative_to(base.resolve()))
        except (OSError, ValueError):
            pass
    return f"[EXTERNAL]/{path.name}"


def find_executable(*names: str) -> Optional[str]:
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


def resolve_iqtree(args: argparse.Namespace) -> str:
    if args.iqtree_executable is not None:
        explicit = args.iqtree_executable.expanduser().resolve()
        if not explicit.is_file() or not os.access(explicit, os.X_OK):
            raise FileNotFoundError(f"--iqtree-executable is not executable: {explicit}")
        return str(explicit)
    tools_dir = Path(__file__).resolve().with_name("tools")
    local_candidates = [tools_dir / "iqtree3"]
    local_candidates.extend(sorted(tools_dir.glob("iqtree-*-macOS/bin/iqtree3"), reverse=True))
    for candidate in local_candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    found = find_executable("iqtree3-omp", "iqtree2-omp", "iqtree-omp", "iqtree3", "iqtree2", "iqtree")
    if found is None:
        raise RuntimeError(
            "IQ-TREE was not found. On macOS run ./install_iqtree_macos.sh, "
            "or provide --iqtree-executable PATH."
        )
    return found


def require_dependencies(args: argparse.Namespace) -> dict[str, str]:
    tools = {
        "mafft": find_executable("mafft"), "trimal": find_executable("trimal"),
        "iqtree": resolve_iqtree(args),
        "blastn": find_executable("blastn"), "makeblastdb": find_executable("makeblastdb"),
        "Rscript": find_executable("Rscript"),
    }
    missing = [name for name, path in tools.items() if path is None]
    if missing:
        formula_for_tool = {
            "mafft": "mafft", "trimal": "trimal", "iqtree": "iqtree3",
            "blastn": "blast", "makeblastdb": "blast", "Rscript": "r",
        }
        formulae = list(dict.fromkeys(formula_for_tool[name] for name in missing))
        raise RuntimeError(
            "Missing external tools: " + ", ".join(missing) + ". On macOS with Homebrew: "
            "brew install " + " ".join(formulae) + "; then install R packages with "
            "Rscript install_R_packages.R."
        )
    if not args.r_script.is_file():
        raise FileNotFoundError(f"R analysis script not found: {args.r_script}")
    check = subprocess.run(
        [tools["Rscript"], "-e", "p<-c('ape','phangorn','phytools');q(status=ifelse(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)),0,2))"],
        check=False,
    )
    if check.returncode:
        raise RuntimeError("Required R packages are missing. Run: Rscript install_R_packages.R")
    return {key: str(value) for key, value in tools.items()}


def run(command: list[str], log: Path, stdout_path: Optional[Path] = None) -> None:
    started = time.perf_counter()
    if stdout_path:
        with stdout_path.open("w", encoding="utf-8", newline="\n") as out:
            result = subprocess.run(command, stdout=out, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(command, capture_output=True, text=True)
    with log.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"$ {shlex.join(command)}\n")
        if not stdout_path:
            handle.write(result.stdout)
        handle.write(result.stderr)
        handle.write(f"\n[elapsed_seconds={time.perf_counter()-started:.3f}; exit={result.returncode}]\n")
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {shlex.join(command)}; see {log}")


def extract_targets(
    summary_path: Path, fasta_path: Path, trait_col: str, out: Path, min_len: int,
    excluded_name_keywords: Optional[list[str]] = None,
    keyword_case_sensitive: bool = False,
) -> tuple[pd.DataFrame, list[tuple[str, str]]]:
    summary = pd.read_csv(summary_path, sep="\t", encoding="utf-8-sig")
    required = {"sh_id", trait_col}
    if not required.issubset(summary.columns):
        raise ValueError(f"Summary lacks columns: {sorted(required - set(summary.columns))}")
    summary[trait_col] = pd.to_numeric(summary[trait_col], errors="coerce")
    summary = summary[summary[trait_col].map(math.isfinite)].drop_duplicates("sh_id")
    wanted = set(summary["sh_id"].astype(str))
    source = {sh_id(h): (h, s) for h, s in read_fasta(fasta_path) if sh_id(h) in wanted}
    excluded_name_keywords = excluded_name_keywords or []
    rows, extracted, kept = [], [], []
    for row in summary.itertuples(index=False):
        sid = str(getattr(row, "sh_id"))
        if sid not in source:
            rows.append({"sh_id": sid, "status": "missing_from_input_fasta", "sequence_length": 0})
            continue
        header, seq = source[sid]
        length = called_length(seq)
        summary_species = getattr(row, "species", None)
        species = (
            str(summary_species).strip().replace(" ", "_")
            if summary_species is not None and pd.notna(summary_species)
            else fasta_species(header)
        )
        count = fasta_n(header)
        # Brackets are Newick comment delimiters, so inference uses a safe label;
        # the exact requested display label is restored only when plotting.
        label = safe_label(f"{sid}_{species}_n{count}")
        display_label = safe_label(f"{sid}_{species}") + f"_[n={count}]"
        haystack = header if keyword_case_sensitive else header.casefold()
        matched_keyword = next((
            keyword for keyword in excluded_name_keywords
            if (keyword if keyword_case_sensitive else keyword.casefold()) in haystack
        ), "")
        if matched_keyword:
            status = "excluded_name_keyword"
        elif length < min_len:
            status = "excluded_below_minimum_length"
        else:
            status = "included"
        summary_genus = getattr(row, "genus", None)
        genus = str(summary_genus).strip() if summary_genus is not None and pd.notna(summary_genus) else "Unknown"
        rows.append({"sh_id": sid, "status": status, "sequence_length": length,
                     "tip_label": label, "display_label": display_label,
                     "species":species.replace("_", " "), "genus":genus,
                     "n_unique_occurrences":count,
                     "matched_exclusion_keyword":matched_keyword})
        extracted.append((label, seq))
        if status == "included":
            kept.append((label, seq))
    report = pd.DataFrame(rows)
    write_fasta(extracted, out / "01_target_sequences_extracted.fasta")
    write_fasta(kept, out / "02_target_sequences_length_filtered.fasta")
    report.to_csv(out / "target_sequence_filter_report.tsv", sep="\t", index=False)
    if len(kept) < 4:
        raise RuntimeError(f"Only {len(kept)} target sequences remain; at least four are required")
    traits = summary.merge(
        report[report.status == "included"], on="sh_id", suffixes=("_summary", "")
    )
    traits = traits[["sh_id", "tip_label", "display_label", "species", "genus",
                     "n_unique_occurrences", trait_col]].rename(columns={trait_col: "trait_value"})
    traits.to_csv(out / "tip_traits.tsv", sep="\t", index=False)
    return traits, kept


def write_extra_trait_table(
    base_traits: pd.DataFrame,
    dprime_path: Path,
    z_column: str,
    output_path: Path,
) -> Path:
    dprime = pd.read_csv(dprime_path, sep="\t", encoding="utf-8-sig")
    required = {"sh_id", z_column}
    if not required.issubset(dprime.columns):
        raise ValueError(f"Trait table lacks columns: {sorted(required - set(dprime.columns))}: {dprime_path}")
    dprime = dprime[["sh_id", z_column]].copy()
    dprime[z_column] = pd.to_numeric(dprime[z_column], errors="coerce")
    merged = base_traits.drop(columns=["trait_value"], errors="ignore").merge(
        dprime, on="sh_id", how="inner"
    )
    merged = merged[merged[z_column].map(math.isfinite)].copy()
    if len(merged) < 4:
        raise RuntimeError(
            f"Only {len(merged)} finite trait values were found in {dprime_path}; at least four are required"
        )
    merged = merged.rename(columns={z_column: "trait_value"})
    merged.to_csv(output_path, sep="\t", index=False)
    return output_path


def unite_stream(path: Path) -> tuple[TextIO, object]:
    if path.name.endswith((".tgz", ".tar.gz")):
        archive = tarfile.open(path, "r:gz")
        members = [m for m in archive.getmembers() if m.isfile() and re.search(r"(?i)\.(fa|fas|fasta)$", m.name)]
        preferred = [m for m in members if "_dev" not in m.name]
        member = (preferred or members)[0]
        raw = archive.extractfile(member)
        if raw is None:
            archive.close(); raise RuntimeError("Could not read UNITE FASTA member")
        return (iter_decode(raw), archive)
    opener = gzip.open if path.suffix == ".gz" else open
    handle = opener(path, "rt", encoding="utf-8", errors="replace")
    return handle, handle


def iter_decode(raw) -> Iterator[str]:
    for line in raw:
        yield line.decode("utf-8", "replace")


def fasta_from_lines(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    header = None; chunks: list[str] = []
    for raw in lines:
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                yield header, clean_seq("".join(chunks))
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line)
    if header is not None:
        yield header, clean_seq("".join(chunks))


def taxonomy(header: str) -> dict[str, str]:
    mapping = {}
    for rank, value in re.findall(r"(?:^|;)([kpcofgs])__([^;|]+)", header):
        mapping[{"k":"kingdom","p":"phylum","c":"class","o":"order","f":"family","g":"genus","s":"species"}[rank]] = value.replace("_", " ")
    parts = header.split("|")
    mapping["accession"] = parts[1] if len(parts) > 1 else "Unknown"
    mapping["sh_id"] = sh_id(header)
    return mapping


def platform_name(header: str) -> str:
    text = header.casefold()
    patterns = [("Sanger", r"sanger|capillary"), ("PacBio", r"pacbio|sequel|revio"),
                ("Illumina", r"illumina|miseq|hiseq|novaseq|nextseq"),
                ("DNBSEQ", r"dnbseq|bgiseq|mgi"), ("Ion Torrent", r"ion.?torrent"),
                ("454 Roche", r"454|roche")]
    for label, pattern in patterns:
        if re.search(pattern, text): return label
    return "Unknown"


def identified_species(value: str) -> bool:
    text = str(value).strip().casefold()
    return bool(text and text not in MISSING_TAXA and not re.search(r"\b(?:sp\.?|cf\.?|aff\.?|uncultured|unidentified|environmental)\b", text))


def identified_genus(value: str) -> bool:
    text = str(value).strip().casefold()
    return bool(text and not re.search(
        r"unknown|unidentified|uncultured|incertae\s+sedis|environmental|\bgen\.?\b|^fungi$",
        text,
    ))


def build_outgroup_candidates(
    unite: Path, target_ids: set[str], target_taxon: str, outgroup_min_len: int,
    max_ambig: float, required_shared_rank: str, minimum_genus_sh_count: int,
    out: Path,
) -> tuple[pd.DataFrame, dict[str, str]]:
    lines, closer = unite_stream(unite)
    target_references, candidates = [], []
    ranks = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
    try:
        for header, seq in fasta_from_lines(lines):
            tax = taxonomy(header); sid = tax["sh_id"]
            target_match = next(
                (rank for rank in ranks if tax.get(rank, "").strip().casefold() == target_taxon.casefold()),
                None,
            )
            if target_match:
                target_references.append(tax)
                continue
            if sid in target_ids:
                continue
            order_name = tax.get("order", "").strip()
            # An outgroup must be positively assigned outside the target order;
            # missing/incertae-sedis order labels are not sufficient evidence.
            if not order_name or "incertae sedis" in order_name.casefold(): continue
            if order_name.casefold() == target_taxon.casefold(): continue
            if not identified_species(tax.get("species", "")): continue
            if not identified_genus(tax.get("genus", "")): continue
            if called_length(seq) < outgroup_min_len or ambiguous_fraction(seq) > max_ambig: continue
            candidates.append((header, seq, tax))
    finally:
        closer.close()
    if not target_references:
        raise RuntimeError(f"Could not infer the UNITE lineage of target taxon {target_taxon!r}")
    if not candidates:
        raise RuntimeError("No species-identified, non-target UNITE outgroup candidates passed filters")
    target_lineage = {}
    for rank in ranks:
        common = Counter(t.get(rank) for t in target_references if t.get(rank)).most_common(1)
        target_lineage[rank] = common[0][0] if common else ""
    required_value = target_lineage.get(required_shared_rank, "")
    if not required_value:
        raise RuntimeError(
            f"Target {required_shared_rank} could not be inferred for {target_taxon!r}; "
            "cannot enforce taxonomic proximity for outgroup selection"
        )
    genus_sh_counts: Counter[str] = Counter()
    genus_to_sh: defaultdict[str, set[str]] = defaultdict(set)
    for _, _, tax in candidates:
        genus = tax.get("genus", "").strip()
        if genus:
            genus_to_sh[genus.casefold()].add(tax["sh_id"])
    genus_sh_counts.update({genus: len(ids) for genus, ids in genus_to_sh.items()})
    records, rows, eligible_records = [], [], []
    for i, (header, seq, tax) in enumerate(candidates, 1):
        cid = f"OGC{i:07d}"
        platform = platform_name(header)
        genus = tax.get("genus", "").strip()
        genus_count = genus_sh_counts.get(genus.casefold(), 0)
        same_required_rank = tax.get(required_shared_rank, "").strip().casefold() == required_value.casefold()
        genus_supported = genus_count >= minimum_genus_sh_count
        eligible = same_required_rank and genus_supported
        tier = 0 if tax.get("class", "").casefold() == target_lineage.get("class", "").casefold() else 1 if tax.get("phylum", "").casefold() == target_lineage.get("phylum", "").casefold() else 2
        rows.append({"candidate_id":cid, **tax, "sequence_length":called_length(seq),
                     "ambiguous_fraction":ambiguous_fraction(seq), "platform":platform,
                     "platform_rank":PLATFORM_ORDER[platform], "phylogenetic_tier":tier,
                     "required_shared_rank":required_shared_rank,
                     "target_rank_value":required_value,
                     "same_required_rank":same_required_rank,
                     "genus_unique_sh_count":genus_count,
                     "minimum_genus_sh_count":minimum_genus_sh_count,
                     "passes_genus_support":genus_supported,
                     "eligible_for_blast":eligible,
                     "original_header":header, "sequence":seq})
        if eligible:
            records.append((cid, seq))
            eligible_records.append(rows[-1])
    report = pd.DataFrame(rows)
    report.drop(columns="sequence").to_csv(out / "outgroup_candidate_filter_report.tsv", sep="\t", index=False)
    if not eligible_records:
        raise RuntimeError(
            f"No outgroup candidate shares target {required_shared_rank}={required_value!r} "
            f"and belongs to a genus represented by at least {minimum_genus_sh_count} unique SHs"
        )
    write_fasta(records, out / "03_unite_outgroup_candidates.fasta")
    lineage_info = {f"target_{rank}": value for rank, value in target_lineage.items()}
    lineage_info.update({
        "required_shared_rank": required_shared_rank,
        "required_shared_rank_value": required_value,
        "minimum_outgroup_genus_sh_count": str(minimum_genus_sh_count),
        "outgroup_minimum_length": str(outgroup_min_len),
        "n_candidates_before_taxonomic_and_genus_filter": str(len(candidates)),
        "n_candidates_eligible_for_blast": str(len(eligible_records)),
    })
    return pd.DataFrame(eligible_records), lineage_info


def choose_outgroup(
    candidates: pd.DataFrame, query: Path, tools: dict[str,str], out: Path,
    jobs: int, minimum_hsp_coverage: float, minimum_median_identity: float,
    maximum_median_identity: float, target_median_identity: float,
    maximum_identity_to_any_target: float, minimum_median_query_coverage: float,
    minimum_query_hit_fraction: float, log: Path,
) -> pd.Series:
    dbprefix = out / "tmp" / "outgroup_db"; dbprefix.parent.mkdir(exist_ok=True)
    run([tools["makeblastdb"], "-in", str(out / "03_unite_outgroup_candidates.fasta"), "-dbtype", "nucl", "-parse_seqids", "-out", str(dbprefix)], log)
    blast = out / "outgroup_blast_hits.tsv"
    fields = "qseqid sseqid pident length qlen slen evalue bitscore qcovs"
    run([tools["blastn"], "-query", str(query), "-db", str(dbprefix), "-outfmt", f"6 {fields}", "-max_target_seqs", "500", "-num_threads", str(jobs)], log, blast)
    hits = pd.read_csv(blast, sep="\t", names=fields.split())
    hits = hits[hits.qcovs >= minimum_hsp_coverage]
    if hits.empty: raise RuntimeError("No outgroup BLAST hit passed minimum query coverage")
    best_hsp = hits.sort_values("bitscore", ascending=False).drop_duplicates(["qseqid","sseqid"])
    stats = best_hsp.groupby("sseqid").agg(
        n_query_hits=("qseqid","nunique"), max_bitscore=("bitscore","max"),
        median_bitscore=("bitscore","median"), median_identity=("pident","median"),
        maximum_identity=("pident","max"), median_query_coverage=("qcovs","median"),
        minimum_evalue=("evalue","min"),
    ).reset_index().rename(columns={"sseqid":"candidate_id"})
    ranked = candidates.drop(columns="sequence").merge(stats, on="candidate_id", how="inner")
    n_target_queries = sum(1 for _ in read_fasta(query))
    ranked["query_hit_fraction"] = ranked.n_query_hits / n_target_queries
    ranked["passes_median_identity_window"] = ranked.median_identity.between(
        minimum_median_identity, maximum_median_identity, inclusive="both"
    )
    ranked["passes_maximum_identity_limit"] = ranked.maximum_identity <= maximum_identity_to_any_target
    ranked["passes_median_query_coverage"] = ranked.median_query_coverage >= minimum_median_query_coverage
    ranked["passes_query_hit_fraction"] = ranked.query_hit_fraction >= minimum_query_hit_fraction
    ranked["goldilocks_eligible"] = (
        ranked.passes_median_identity_window
        & ranked.passes_maximum_identity_limit
        & ranked.passes_median_query_coverage
        & ranked.passes_query_hit_fraction
    )
    ranked["distance_from_target_median_identity"] = (
        ranked.median_identity - target_median_identity
    ).abs()
    pool = ranked[ranked.goldilocks_eligible].copy()
    if pool.empty:
        ranked.to_csv(out / "outgroup_candidate_ranking.tsv", sep="\t", index=False)
        raise RuntimeError(
            "No same-class outgroup candidate passed the Goldilocks distance/coverage "
            "criteria. Inspect outgroup_candidate_ranking.tsv and adjust the explicit "
            "--outgroup-* thresholds if scientifically justified."
        )
    pool = pool.sort_values(
        ["distance_from_target_median_identity", "platform_rank",
         "median_query_coverage", "query_hit_fraction", "genus_unique_sh_count",
         "ambiguous_fraction", "sequence_length"],
        ascending=[True, True, False, False, False, True, False],
    )
    selected = pool.iloc[0]
    ranked["selected"] = ranked.candidate_id.eq(selected.candidate_id)
    ranked = ranked.sort_values(
        ["selected", "goldilocks_eligible", "distance_from_target_median_identity",
         "platform_rank", "median_query_coverage"],
        ascending=[False, False, True, True, False],
    )
    ranked.to_csv(out / "outgroup_candidate_ranking.tsv", sep="\t", index=False)
    chosen = candidates.set_index("candidate_id").loc[selected.candidate_id].copy()
    for key, value in selected.items():
        if key != "candidate_id":
            chosen[key] = value
    return chosen


def aligned_pair_identity(seq_a: str, seq_b: str) -> tuple[float, int]:
    comparable = [(a, b) for a, b in zip(seq_a, seq_b) if a in "ACGT" and b in "ACGT"]
    overlap = len(comparable)
    if overlap == 0:
        return math.nan, 0
    matches = sum(a == b for a, b in comparable)
    return 100.0 * matches / overlap, overlap


def screen_targets_by_similarity(
    aligned_path: Path, targets: list[tuple[str, str]], outgroup_label: str,
    out: Path, minimum_overlap: int, minimum_comparisons: int,
    identity_margin: float, mad_multiplier: float,
) -> tuple[list[tuple[str, str]], pd.DataFrame]:
    aligned = dict(read_fasta(aligned_path))
    target_labels = [label for label, _ in targets]
    missing = [label for label in [*target_labels, outgroup_label] if label not in aligned]
    if missing:
        raise RuntimeError(f"Pre-screen alignment is missing {len(missing)} expected sequence(s)")

    outgroup_comparisons = []
    identity_to_outgroup: dict[str, float] = {}
    overlap_to_outgroup: dict[str, int] = {}
    for label in target_labels:
        identity, overlap = aligned_pair_identity(aligned[outgroup_label], aligned[label])
        identity_to_outgroup[label] = identity
        overlap_to_outgroup[label] = overlap
        if overlap >= minimum_overlap and math.isfinite(identity):
            outgroup_comparisons.append(identity)
    if len(outgroup_comparisons) < minimum_comparisons:
        raise RuntimeError(
            "Too few outgroup-to-target comparisons passed the overlap threshold for "
            "target taxonomic screening"
        )
    outgroup_mean = sum(outgroup_comparisons) / len(outgroup_comparisons)
    outgroup_baseline = outgroup_mean + identity_margin
    rows = []
    for label in target_labels:
        target_comparisons = []
        overlaps = []
        for other in target_labels:
            if other == label:
                continue
            identity, overlap = aligned_pair_identity(aligned[label], aligned[other])
            if overlap >= minimum_overlap and math.isfinite(identity):
                target_comparisons.append(identity)
                overlaps.append(overlap)
        mean_identity = (
            sum(target_comparisons) / len(target_comparisons)
            if target_comparisons else math.nan
        )
        rows.append({
            "tip_label": label,
            "sh_id": sh_id(label),
            "mean_identity_to_other_targets_percent": mean_identity,
            "n_target_comparisons": len(target_comparisons),
            "median_pairwise_overlap_bp": float(pd.Series(overlaps).median()) if overlaps else math.nan,
            "identity_to_outgroup_percent": identity_to_outgroup[label],
            "overlap_with_outgroup_bp": overlap_to_outgroup[label],
            "outgroup_mean_identity_to_targets_percent": outgroup_mean,
            "outgroup_baseline_with_margin_percent": outgroup_baseline,
            "identity_margin_percentage_points": identity_margin,
            "minimum_pairwise_overlap_bp": minimum_overlap,
            "minimum_comparisons": minimum_comparisons,
            "outgroup_tip_label": outgroup_label,
        })
    report = pd.DataFrame(rows)
    eligible_means = report.loc[
        report.n_target_comparisons >= minimum_comparisons,
        "mean_identity_to_other_targets_percent",
    ].dropna()
    target_median = float(eligible_means.median())
    target_mad = float((eligible_means - target_median).abs().median())
    robust_low_outlier_cutoff = target_median - mad_multiplier * target_mad
    effective_cutoff = min(outgroup_baseline, robust_low_outlier_cutoff)
    report["target_mean_identity_median_percent"] = target_median
    report["target_mean_identity_mad_percent"] = target_mad
    report["mad_multiplier"] = mad_multiplier
    report["robust_low_outlier_cutoff_percent"] = robust_low_outlier_cutoff
    report["effective_exclusion_cutoff_percent"] = effective_cutoff
    sufficient = report.n_target_comparisons >= minimum_comparisons
    report["retained"] = ~sufficient | (
        report.mean_identity_to_other_targets_percent > effective_cutoff
    )
    report["status"] = "retained_above_conservative_dual_cutoff"
    report.loc[~sufficient, "status"] = "retained_insufficient_comparisons"
    report.loc[sufficient & ~report.retained, "status"] = (
        "excluded_below_outgroup_baseline_and_robust_low_outlier_cutoff"
    )
    kept_labels = set(report.loc[report.retained, "tip_label"])
    report = report.sort_values(
        ["retained", "mean_identity_to_other_targets_percent", "sh_id"],
        ascending=[True, True, True], na_position="last",
    )
    report.to_csv(out / "target_taxonomic_similarity_screen.tsv", sep="\t", index=False)
    summary = pd.DataFrame([{
        "n_targets_before_screen": len(targets),
        "n_targets_retained": len(kept_labels),
        "n_targets_excluded": len(targets) - len(kept_labels),
        "outgroup_mean_identity_to_targets_percent": outgroup_mean,
        "identity_margin_percentage_points": identity_margin,
        "outgroup_baseline_with_margin_percent": outgroup_baseline,
        "target_mean_identity_median_percent": target_median,
        "target_mean_identity_mad_percent": target_mad,
        "mad_multiplier": mad_multiplier,
        "robust_low_outlier_cutoff_percent": robust_low_outlier_cutoff,
        "effective_exclusion_cutoff_percent": effective_cutoff,
        "minimum_pairwise_overlap_bp": minimum_overlap,
        "minimum_comparisons": minimum_comparisons,
        "criterion": "exclude only when target mean identity is at or below both the outgroup baseline plus margin and the target-median minus MAD-multiplier cutoff; retain if comparisons are insufficient",
    }])
    summary.to_csv(out / "target_taxonomic_similarity_screen_summary.tsv", sep="\t", index=False)
    kept = [(label, seq) for label, seq in targets if label in kept_labels]
    if len(kept) < 4:
        raise RuntimeError(f"Taxonomic similarity screening retained only {len(kept)} targets")
    write_fasta(kept, out / "02b_target_sequences_taxonomically_screened.fasta")
    return kept, report


def fasta_statistics(path: Path) -> dict[str, object]:
    records = list(read_fasta(path))
    lengths = [called_length(seq) for _, seq in records]
    raw_lengths = [len(seq) for _, seq in records]
    aligned = bool(records) and len(set(raw_lengths)) == 1
    cells = sum(raw_lengths)
    gaps = sum(seq.count("-") + seq.count(".") for _, seq in records)
    ambiguous = sum(sum(base not in "ACGT-?." for base in seq) for _, seq in records)
    return {
        "n_sequences": len(records),
        "total_called_bases": sum(lengths),
        "minimum_called_length": min(lengths) if lengths else 0,
        "median_called_length": float(pd.Series(lengths).median()) if lengths else 0,
        "maximum_called_length": max(lengths) if lengths else 0,
        "is_equal_length_alignment": aligned,
        "alignment_columns": raw_lengths[0] if aligned else "",
        "total_matrix_cells": cells if aligned else "",
        "gap_cells": gaps if aligned else "",
        "ambiguous_cells": ambiguous if aligned else "",
        "gap_fraction": gaps / cells if aligned and cells else "",
        "ambiguous_fraction": ambiguous / cells if aligned and cells else "",
    }


def write_alignment_column_occupancy(path: Path, output: Path) -> pd.DataFrame:
    records = list(read_fasta(path))
    if not records:
        raise RuntimeError(f"No sequences found in alignment: {path}")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise RuntimeError(f"Sequences have unequal aligned lengths: {path}")
    n_sequences = len(records)
    rows = []
    for column_index, column in enumerate(zip(*(seq for _, seq in records)), start=1):
        n_gap = sum(base in "-." for base in column)
        n_called = sum(base.upper() in "ACGT" for base in column)
        n_non_gap = n_sequences - n_gap
        rows.append({
            "alignment_column": column_index,
            "n_sequences": n_sequences,
            "n_non_gap_sequences": n_non_gap,
            "non_gap_fraction": n_non_gap / n_sequences,
            "n_called_acgt_sequences": n_called,
            "called_acgt_fraction": n_called / n_sequences,
            "n_gap_sequences": n_gap,
        })
    result = pd.DataFrame(rows)
    result.to_csv(output, sep="\t", index=False)
    return result


def trim_alignment_terminals(
    input_path: Path,
    output_path: Path,
    report_path: Path,
    minimum_occupancy: float,
) -> dict[str, object]:
    records = list(read_fasta(input_path))
    if not records:
        raise RuntimeError(f"No sequences found in alignment: {input_path}")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise RuntimeError(f"Sequences have unequal aligned lengths: {input_path}")
    n_sequences = len(records)
    input_columns = lengths.pop()
    required_sequences = math.ceil(
        minimum_occupancy * n_sequences - 1e-12
    )
    called_counts = [
        sum(base.upper() in "ACGT" for base in column)
        for column in zip(*(seq for _, seq in records))
    ]
    eligible = [
        index for index, count in enumerate(called_counts)
        if count >= required_sequences
    ]
    if not eligible:
        raise RuntimeError(
            "No MAFFT alignment column meets "
            f"--terminal-trim-minimum-occupancy {minimum_occupancy} "
            f"({required_sequences}/{n_sequences} sequences with A/C/G/T). "
            "Inspect the alignment or lower the option explicitly."
        )
    left = eligible[0]
    right = eligible[-1]
    trimmed_records = [
        (label, sequence[left:right + 1]) for label, sequence in records
    ]
    write_fasta(trimmed_records, output_path)
    summary = {
        "input_file": input_path.name,
        "output_file": output_path.name,
        "n_sequences": n_sequences,
        "minimum_occupancy_fraction": minimum_occupancy,
        "minimum_sequences_required": required_sequences,
        "input_alignment_columns": input_columns,
        "left_boundary_input_column_1based": left + 1,
        "right_boundary_input_column_1based": right + 1,
        "columns_removed_from_left": left,
        "columns_removed_from_right": input_columns - right - 1,
        "output_alignment_columns": right - left + 1,
        "left_boundary_called_sequences": called_counts[left],
        "right_boundary_called_sequences": called_counts[right],
        "boundary_character_definition": "unambiguous_A_C_G_T",
        "internal_columns_preserved_before_trimAl": True,
    }
    pd.DataFrame([summary]).to_csv(report_path, sep="\t", index=False)
    return summary


def iqtree_report_statistics(text: str) -> dict[str, object]:
    patterns = {
        "selected_substitution_model": (
            r"Best-fit model according to BIC:\s*(\S+)",
            r"Best-fit model:\s*(\S+)\s+chosen according to BIC",
        ),
        "parsimony_informative_sites": (
            r"Number of parsimony informative sites:\s*(\d+)",
            r"(\d+) parsimony-informative",
        ),
        "singleton_sites": (r"(\d+) singleton sites",),
        "constant_sites": (r"Number of constant sites:\s*(\d+)", r"(\d+) constant sites"),
        "distinct_site_patterns": (r"Number of distinct site patterns:\s*(\d+)", r"(\d+) distinct patterns"),
        "optimal_log_likelihood": (r"Log-likelihood of the tree:\s*([-+0-9.eE]+)", r"Optimal log-likelihood:\s*([-+0-9.eE]+)"),
        "total_tree_length": (r"Total tree length \(sum of branch lengths\):\s*([-+0-9.eE]+)", r"Total tree length:\s*([-+0-9.eE]+)"),
    }
    result = {}
    for metric, alternatives in patterns.items():
        match = next((m for pattern in alternatives if (m := re.search(pattern, text))), None)
        result[metric] = match.group(1) if match else "not_reported"
    return result


def write_publication_summary(
    out: Path, stages: list[tuple[str, Path]], iqtext: str, chosen: pd.Series,
    lineage: dict[str, str], bootstrap: int, iqtree_jobs: int,
    trimal_gap_threshold: float, trimal_minimum_retained_percent: float,
) -> None:
    rows = []
    for stage, path in stages:
        if not path.is_file():
            continue
        for metric, value in fasta_statistics(path).items():
            rows.append({"category":"sequence_stage", "item":stage, "metric":metric,
                         "value":value, "unit":"", "source_file":path.name})
    for metric, value in iqtree_report_statistics(iqtext).items():
        rows.append({"category":"maximum_likelihood", "item":"IQ-TREE", "metric":metric,
                     "value":value, "unit":"", "source_file":"ml_iqtree.iqtree"})
    for metric, value in {
        "standard_nonparametric_bootstrap_replicates": bootstrap,
        "iqtree_threads_used": iqtree_jobs,
    }.items():
        rows.append({"category":"maximum_likelihood", "item":"IQ-TREE", "metric":metric,
                     "value":value, "unit":"", "source_file":"run_configuration.tsv"})
    for metric, value, unit in (
        ("gap_score_threshold", trimal_gap_threshold, "fraction_non_gap_required"),
        ("minimum_alignment_retained", trimal_minimum_retained_percent, "percent"),
    ):
        rows.append({"category":"alignment_trimming", "item":"trimAl_manual_gap_preserving",
                     "metric":metric, "value":value, "unit":unit,
                     "source_file":"run_configuration.tsv"})
    for metric, value in lineage.items():
        rows.append({"category":"outgroup_selection", "item":"target_lineage", "metric":metric,
                     "value":value, "unit":"", "source_file":"UNITE"})
    for metric in ("sh_id", "species", "genus", "family", "order", "class", "phylum",
                   "sequence_length", "ambiguous_fraction", "platform", "genus_unique_sh_count",
                   "median_identity", "maximum_identity", "median_query_coverage",
                   "query_hit_fraction", "distance_from_target_median_identity",
                   "max_bitscore", "goldilocks_eligible"):
        if metric in chosen.index:
            rows.append({"category":"outgroup_selection", "item":"selected_outgroup", "metric":metric,
                         "value":chosen[metric], "unit":"", "source_file":"selected_outgroup.tsv"})
    screen_summary_path = out / "target_taxonomic_similarity_screen_summary.tsv"
    if screen_summary_path.is_file():
        screen_summary = pd.read_csv(screen_summary_path, sep="\t").iloc[0]
        for metric, value in screen_summary.items():
            rows.append({"category":"target_taxonomic_screen",
                         "item":"mean_pairwise_identity_vs_outgroup_baseline",
                         "metric":metric, "value":value, "unit":"",
                         "source_file":screen_summary_path.name})
    pd.DataFrame(rows).to_csv(out / "publication_analysis_summary.tsv", sep="\t", index=False)


def main() -> int:
    args = parse_args(); started = RUN_STARTED_MONOTONIC; started_at = RUN_STARTED_AT
    project_dir = Path(__file__).resolve().parent
    args.output_dir = args.output_root / args.specificity_run_name
    args.data_dir = args.data_dir.resolve(); args.output_root = args.output_root.resolve(); args.output_dir = args.output_dir.resolve(); args.output_dir.mkdir(parents=True, exist_ok=True)
    args.r_script = args.r_script.resolve()
    if args.iqtree_executable is not None:
        args.iqtree_executable = args.iqtree_executable.expanduser().resolve()
    jobs = args.jobs or (os.cpu_count() or 1)
    specificity_inputs = resolve_specificity_inputs(args, project_dir)
    summary = (args.summary or specificity_inputs.get("summary") or one_match(
        args.data_dir, "*selected_fungi_integrated_summary.tsv", "integrated summary"
    )).resolve()
    two_dp_table = (args.two_dp_table or specificity_inputs.get("two_dp_table") or one_match(
        args.data_dir, "*sampling_unit_label_shuffle_2dp_z.tsv", "2DP z-score table"
    )).resolve()
    input_fasta = (args.input_fasta or specificity_inputs.get("input_fasta") or one_match(
        args.data_dir, "*root_host_sequences_its1_and_its2.fasta", "ITS1+ITS2 FASTA"
    )).resolve()
    continent_dprime = None
    continent_two_dp = None
    unite = (args.unite or one_match(args.data_dir, "sh_general_release*.tgz", "UNITE archive")).resolve()
    log = args.output_dir / "workflow.log"
    command = shlex.join([sys.executable, *sys.argv])
    log.write_text(f"Started: {started_at.isoformat()}\nCommand: {command}\n", encoding="utf-8")
    (args.output_dir / "run_command.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + command + "\n", encoding="utf-8")
    config = vars(args) | {
        "summary_resolved":summary, "two_dp_table_resolved":two_dp_table,
        "input_fasta_resolved":input_fasta, "unite_resolved":unite,
        "continent_dprime_resolved": continent_dprime or "",
        "continent_two_dp_resolved": continent_two_dp or "",
        **{f"{key}_resolved": value for key, value in specificity_inputs.items()},
        "jobs_resolved":jobs, "sequence_type":"non-protein-coding ribosomal/ITS DNA",
        "ml_bootstrap_type":"standard_nonparametric",
        "trimal_policy":"manual_gap_preserving",
    }
    tools = require_dependencies(args)
    iqtree_version = command_version(tools["iqtree"])
    iqtree_parallel_capable = "single-core" not in iqtree_version.casefold()
    if not iqtree_parallel_capable:
        iqtree_thread_args = ["-nt", "1"]
        iqtree_threads_requested: object = 1
    elif args.iqtree_thread_mode == "auto":
        iqtree_thread_args = ["-nt", "AUTO", "-ntmax", str(jobs)]
        iqtree_threads_requested = f"AUTO (maximum {jobs})"
    else:
        iqtree_thread_args = ["-nt", str(jobs)]
        iqtree_threads_requested = jobs
    config.update({
        "iqtree_executable_resolved": tools["iqtree"],
        "iqtree_version_detected": iqtree_version,
        "iqtree_parallel_capable": iqtree_parallel_capable,
        "iqtree_threads_requested": iqtree_threads_requested,
        "iqtree_threads_maximum": jobs,
    })
    path_settings = {
        "data_dir", "summary", "two_dp_table", "input_fasta", "unite", "output_root", "output_dir",
        "specificity_output_root",
        "r_script", "iqtree_executable", "summary_resolved", "two_dp_table_resolved",
        "input_fasta_resolved", "unite_resolved", "iqtree_executable_resolved",
        "continent_dprime_resolved", "continent_two_dp_resolved",
        "specificity_run_dir_resolved", "specificity_result_dir_resolved",
    }
    config_rows = []
    for key, value in config.items():
        if key in path_settings and value not in (None, ""):
            value = portable_report_path(value, project_dir)
        config_rows.append({"setting":key, "value":value})
    pd.DataFrame(config_rows).to_csv(args.output_dir/"run_configuration.tsv",sep="\t",index=False)
    traits, targets = extract_targets(
        summary, input_fasta, args.trait_column, args.output_dir, args.minimum_length,
        args.exclude_sequence_name_keyword,
        args.sequence_name_keyword_case_sensitive,
    )
    candidates, target_lineage = build_outgroup_candidates(
        unite, set(traits.sh_id), args.target_taxon, args.outgroup_minimum_length,
        args.maximum_ambiguous_fraction, args.outgroup_required_shared_rank,
        args.minimum_outgroup_genus_sh_count, args.output_dir,
    )
    chosen = choose_outgroup(
        candidates, args.output_dir/"02_target_sequences_length_filtered.fasta",
        tools, args.output_dir, jobs, args.minimum_outgroup_query_coverage,
        args.outgroup_minimum_median_identity,
        args.outgroup_maximum_median_identity,
        args.outgroup_target_median_identity,
        args.outgroup_maximum_identity_to_any_target,
        args.outgroup_minimum_median_query_coverage,
        args.outgroup_minimum_query_hit_fraction,
        log,
    )
    outgroup_label = safe_label(f"OUTGROUP_{chosen['sh_id']}_{chosen['species']}")
    write_fasta([(outgroup_label, str(chosen["sequence"]))], args.output_dir/"03b_selected_outgroup.fasta")
    pd.DataFrame([{**chosen.drop(labels="sequence").to_dict(), "tip_label":outgroup_label, "selection_notes":"outside target taxon; same class; genus supported by the required unique-SH count; species identified; length/ambiguity passed; broad query coverage; median identity in the configured Goldilocks window; no excessively close target hit; closest to target median identity; platform and sequence quality used as tie-breakers"}]).to_csv(args.output_dir/"selected_outgroup.tsv",sep="\t",index=False)
    mafft_flags={"linsi":["--localpair","--maxiterate","1000"],"einsi":["--genafpair","--maxiterate","1000"],"ginsi":["--globalpair","--maxiterate","1000"],"auto":["--auto"]}[args.mafft_algorithm]
    preliminary_unaligned = args.output_dir/"04a_targets_plus_outgroup_pre_screen_unaligned.fasta"
    preliminary_aligned = args.output_dir/"04b_mafft_pre_screen_aligned.fasta"
    write_fasta(targets + [(outgroup_label, str(chosen["sequence"]))], preliminary_unaligned)
    run([tools["mafft"], *mafft_flags, "--thread", str(jobs), str(preliminary_unaligned)], log, preliminary_aligned)
    traits.to_csv(args.output_dir/"tip_traits_before_taxonomic_screen.tsv", sep="\t", index=False)
    targets, similarity_report = screen_targets_by_similarity(
        preliminary_aligned, targets, outgroup_label, args.output_dir,
        args.misidentification_minimum_pairwise_overlap,
        args.misidentification_minimum_comparisons,
        args.target_vs_outgroup_mean_identity_margin,
        args.misidentification_mad_multiplier,
    )
    retained_labels = {label for label, _ in targets}
    traits = traits[traits.tip_label.isin(retained_labels)].copy()
    traits.to_csv(args.output_dir/"tip_traits.tsv", sep="\t", index=False)
    extra_r_args: list[str] = []
    if continent_dprime is not None and continent_two_dp is not None:
        continent_traits = write_extra_trait_table(
            traits,
            continent_dprime.resolve(),
            "z_standardized_dprime",
            args.output_dir / "tip_traits_continent.tsv",
        )
        extra_r_args.extend([
            str(continent_traits),
            str(continent_two_dp.resolve()),
            "continent",
            "continent-specificity z-d'",
            "z-standardized d'",
        ])
    with_outgroup = targets + [(outgroup_label, str(chosen["sequence"]))]
    write_fasta(with_outgroup, args.output_dir/"04_targets_plus_outgroup_unaligned.fasta")
    aligned=args.output_dir/"05_mafft_aligned.fasta"
    run([tools["mafft"],*mafft_flags,"--thread",str(jobs),str(args.output_dir/"04_targets_plus_outgroup_unaligned.fasta")],log,aligned)
    terminal_trimmed = args.output_dir/"05b_terminal_trimmed_alignment.fasta"
    terminal_trim_summary = trim_alignment_terminals(
        aligned,
        terminal_trimmed,
        args.output_dir/"terminal_trimming_summary.tsv",
        args.terminal_trim_minimum_occupancy,
    )
    with log.open("a", encoding="utf-8") as handle:
        handle.write(
            "\nTerminal alignment trimming: "
            f"columns {terminal_trim_summary['left_boundary_input_column_1based']}-"
            f"{terminal_trim_summary['right_boundary_input_column_1based']} retained; "
            f"{terminal_trim_summary['output_alignment_columns']} of "
            f"{terminal_trim_summary['input_alignment_columns']} columns remain; "
            f"boundary occupancy={terminal_trim_summary['minimum_sequences_required']}/"
            f"{terminal_trim_summary['n_sequences']} sequences.\n"
        )
    trimmed=args.output_dir/"06_trimal_trimmed.fasta"
    trimal_command = [
        tools["trimal"], "-in", str(terminal_trimmed), "-out", str(trimmed),
        "-gt", str(args.trimal_gap_threshold),
        "-htmlout", str(args.output_dir/"trimal_report.html"),
    ]
    if args.trimal_minimum_retained_percent > 0:
        trimal_command.extend(["-cons", str(args.trimal_minimum_retained_percent)])
    run(trimal_command, log)
    occupancy = write_alignment_column_occupancy(
        trimmed, args.output_dir/"trimal_column_occupancy.tsv"
    )
    singleton_columns = int((occupancy["n_non_gap_sequences"] == 1).sum())
    configured_failures = int(
        (occupancy["non_gap_fraction"] + 1e-12 < args.trimal_gap_threshold).sum()
    )
    with log.open("a", encoding="utf-8") as handle:
        handle.write(
            "\ntrimAl occupancy audit: "
            f"{len(occupancy)} retained columns; "
            f"minimum non-gap sequences={int(occupancy['n_non_gap_sequences'].min())}; "
            f"singleton columns={singleton_columns}; "
            f"columns below -gt threshold={configured_failures}.\n"
        )
        if args.trimal_minimum_retained_percent > 0 and configured_failures:
            handle.write(
                "WARNING: --trimal-minimum-retained-percent caused trimAl to retain "
                "columns below --trimal-gap-threshold. Set it to 0 to enforce the "
                "gap threshold without a minimum-retention override.\n"
            )
    mlprefix=args.output_dir/"ml_iqtree"; run([
        tools["iqtree"], "-s", str(trimmed), "-st", "DNA", "-m", "MFP",
        "-b", str(args.bootstrap), *iqtree_thread_args, "-seed", str(args.seed),
        "-o", outgroup_label, "-pre", str(mlprefix), "-redo",
    ], log)
    ml_tree=Path(str(mlprefix)+".contree");
    if not ml_tree.exists(): ml_tree=Path(str(mlprefix)+".treefile")
    iqtext=Path(str(mlprefix)+".iqtree").read_text(encoding="utf-8",errors="replace")
    iqlogtext=Path(str(mlprefix)+".log").read_text(encoding="utf-8",errors="replace")
    thread_matches = re.findall(r"BEST NUMBER OF THREADS:\s*(\d+)", iqlogtext)
    if not thread_matches:
        thread_matches = re.findall(r"Kernel:.*?-\s*(\d+) threads", iqlogtext)
    iqtree_threads_actual = int(thread_matches[-1]) if thread_matches else "not_reported"
    config["iqtree_threads_resolved"] = iqtree_threads_actual
    config_rows = []
    for key, value in config.items():
        if key in path_settings and value not in (None, ""):
            value = portable_report_path(value, project_dir)
        config_rows.append({"setting":key, "value":value})
    pd.DataFrame(config_rows).to_csv(args.output_dir/"run_configuration.tsv",sep="\t",index=False)
    model = (
        re.search(r"Best-fit model according to BIC:\s*(\S+)", iqtext)
        or re.search(r"Best-fit model:\s*(\S+)\s+chosen according to BIC", iqtext)
    )
    pd.DataFrame([{"sequence_type":"non-protein-coding DNA","criterion":"BIC","selected_model":model.group(1) if model else "see_iqtree_report","iqtree_report":Path(str(mlprefix)+".iqtree").name}]).to_csv(args.output_dir/"selected_substitution_model.tsv",sep="\t",index=False)
    write_publication_summary(
        args.output_dir,
        [
            ("target_sequences_extracted", args.output_dir/"01_target_sequences_extracted.fasta"),
            ("target_sequences_length_filtered", args.output_dir/"02_target_sequences_length_filtered.fasta"),
            ("target_sequences_taxonomically_screened", args.output_dir/"02b_target_sequences_taxonomically_screened.fasta"),
            ("eligible_outgroup_candidates", args.output_dir/"03_unite_outgroup_candidates.fasta"),
            ("selected_outgroup", args.output_dir/"03b_selected_outgroup.fasta"),
            ("targets_plus_outgroup_pre_screen_unaligned", preliminary_unaligned),
            ("mafft_pre_screen_alignment", preliminary_aligned),
            ("targets_plus_outgroup_unaligned", args.output_dir/"04_targets_plus_outgroup_unaligned.fasta"),
            ("mafft_alignment", aligned),
            ("terminal_trimmed_alignment", terminal_trimmed),
            ("trimmed_alignment", trimmed),
        ],
        iqtext, chosen, target_lineage, args.bootstrap, iqtree_threads_actual,
        args.trimal_gap_threshold, args.trimal_minimum_retained_percent,
    )
    run([
        tools["Rscript"], str(args.r_script), str(trimmed), str(ml_tree),
        str(args.output_dir/"tip_traits.tsv"), str(two_dp_table), str(args.output_dir),
        outgroup_label, str(args.bootstrap), str(args.signal_randomizations),
        str(args.seed), args.nj_distance_model, str(args.figure_font_scale),
        *extra_r_args,
    ], log)
    versions=[]
    for name,path in tools.items(): versions.append({"category":"external_tool","name":name,"version":command_version(path),"path":portable_report_path(path, project_dir)})
    versions += [{"category":"python","name":"Python","version":platform.python_version(),"path":f"[EXTERNAL]/{Path(sys.executable).name}"}]
    for script in (Path(__file__).resolve(), args.r_script.resolve()):
        versions.append({"category":"script","name":script.name,"version":hashlib.sha256(script.read_bytes()).hexdigest(),"path":portable_report_path(script, project_dir)})
    for dist in importlib_metadata.distributions(): versions.append({"category":"python_package","name":dist.metadata.get("Name","Unknown"),"version":dist.version,"path":"[PYTHON_ENV]"})
    pd.DataFrame(versions).to_csv(args.output_dir/"runtime_versions.tsv",sep="\t",index=False)
    ended=datetime.now(timezone.utc).astimezone(); elapsed=time.perf_counter()-started
    pd.DataFrame([{"start_time":started_at.isoformat(),"end_time":ended.isoformat(),"elapsed_seconds":elapsed,"status":"completed","command":command}]).to_csv(args.output_dir/"run_timing.tsv",sep="\t",index=False)
    print(f"Completed phylogeny workflow: {args.output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except BaseException:
        def cli_value(flag: str, default: Optional[str] = None) -> Optional[str]:
            if flag in sys.argv:
                index = sys.argv.index(flag)
                if index + 1 < len(sys.argv):
                    return sys.argv[index + 1]
            return default
        output_root = Path(cli_value("--output-root", "outputs") or "outputs")
        output = output_root / (cli_value("--specificity-run-name", "failed_run") or "failed_run")
        output = output.resolve()
        if output.exists():
            with (output / "workflow.log").open("a", encoding="utf-8") as handle:
                handle.write("\nFAILED\n" + traceback.format_exc() + "\n")
            pd.DataFrame([{
                "start_time": RUN_STARTED_AT.isoformat(),
                "end_time": datetime.now(timezone.utc).astimezone().isoformat(),
                "elapsed_seconds": time.perf_counter() - RUN_STARTED_MONOTONIC,
                "status": "failed", "command": shlex.join([sys.executable, *sys.argv]),
            }]).to_csv(output / "run_timing.tsv", sep="\t", index=False)
        raise
