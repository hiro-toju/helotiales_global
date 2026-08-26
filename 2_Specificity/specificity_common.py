#!/usr/bin/env python3
"""Host-plant and continent specificity using sampling-unit label shuffling.

Unique-occurrence TSV files are read directly. Site-host sampling units are
reconstructed, and their plant labels are shuffled without replacement across
all sampling units, without spatial blocks. A fungal-by-plant matrix is rebuilt before
global d' and 2DP are calculated. The same block-free procedure independently
shuffles continent labels to evaluate fungus-to-continent specificity.
No matrix-level null algorithm is used.
"""

from __future__ import annotations

import argparse
import csv
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
import tempfile
import time
import urllib.parse
import urllib.request
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Iterable, Iterator, Optional, Union

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "specificity_matplotlib_cache")
)

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.text import Text
import numpy as np
import pandas as pd


UNKNOWN_LABELS = {"", "unknown", "unidentified", "unresolved", "na", "nan", "none"}
TAXONOMIC_RANKS = ("genus", "family", "order", "class", "phylum", "kingdom")
CONFIDENCE_LABELS = ("medium_or_high",)
GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
EXCEL_CELL_LIMIT = 32_767
RUN_STATE: dict[str, object] = {}
FIGURE_FONT_SCALE = 1.0


def set_figure_font_scale(scale: float) -> None:
    """Set the multiplier applied to every text element before figure export."""
    global FIGURE_FONT_SCALE
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("--figure-font-scale must be a finite number greater than 0.")
    FIGURE_FONT_SCALE = float(scale)


def apply_figure_font_scale(fig: plt.Figure) -> None:
    """Apply the configured multiplier once to every text artist in a figure."""
    if not getattr(fig, "_specificity_font_scale_applied", False):
        for text in fig.findobj(match=Text):
            text.set_fontsize(text.get_fontsize() * FIGURE_FONT_SCALE)
        fig._specificity_font_scale_applied = True


def save_figure(fig: plt.Figure, output_path: Path, dpi: int = 300) -> None:
    """Save each figure as a vector PDF."""
    apply_figure_font_scale(fig)
    vector_path = output_path.with_suffix(".pdf")
    fig.savefig(vector_path, format="pdf", bbox_inches="tight")


def parse_args(analysis_level: str) -> argparse.Namespace:
    if analysis_level not in {"sh", "genus"}:
        raise ValueError("analysis_level must be 'sh' or 'genus'")
    default_glob = "*root_occurrences_with_hosts.tsv"
    parser = argparse.ArgumentParser(
        description=(
            f"Host and continent specificity for fungal {analysis_level.upper()} rows, "
            "using global sampling-unit label shuffling without blocks."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument(
        "--globalfungi-results", type=Path, default=None,
        help=(
            "Directory produced by 1_GlobalFungi, e.g. ../1_GlobalFungi/results_helotiales. "
            "The medium-or-higher occurrence TSV and matching FASTA files are "
            "copied into this run's data/ folder before analysis."
        ),
    )
    parser.add_argument("--input-glob", default=default_glob)
    parser.add_argument(
        "--plant-rank", choices=TAXONOMIC_RANKS, default="family",
        help="Plant rank used to aggregate columns (default: family).",
    )
    parser.add_argument(
        "--fungi-selection", choices=("rank", "count"), default="rank",
        help="Select fungi by sampling-unit rank or minimum sampling-unit count (default: rank).",
    )
    parser.add_argument("--fungi-rank-max", type=int, default=50)
    parser.add_argument("--fungi-min-count", type=int, default=1)
    parser.add_argument(
        "--plant-selection", choices=("rank", "count"), default="rank",
        help="Select plant lineages by sampling-unit rank or minimum sampling-unit count (default: rank).",
    )
    parser.add_argument("--plant-rank-max", type=int, default=30)
    parser.add_argument("--plant-min-count", type=int, default=1)
    parser.add_argument("--n-randomizations", type=int, default=100)
    parser.add_argument(
        "--figure-font-scale", type=float, default=1.5,
        help="Multiplier for every figure text size (default: 1.5).",
    )
    parser.add_argument(
        "--jobs", type=int, default=0,
        help=(
            "Worker budget. 0 uses all available CPU cores. "
            "No artificial cap is applied outside IQ-TREE."
        ),
    )
    parser.add_argument("--seed", type=int, default=20260704)
    parser.add_argument(
        "--r-script", type=Path,
        default=Path(__file__).with_name("specificity_sample_unit_label_shuffle.R"),
        help="R helper that globally shuffles sampling-unit labels without blocks.",
    )
    parser.add_argument(
        "--world-geojson", type=Path, default=None,
        help=(
            "Natural Earth GeoJSON for SH maps "
            "(default: assets/ne_110m_admin_0_countries.geojson)."
        ),
    )
    parser.add_argument(
        "--taxonomy", type=Path, default=None,
        help="Cached genus taxonomy TSV (default: DATA_DIR/plant_taxonomy.tsv).",
    )
    parser.add_argument(
        "--taxonomy-overrides", type=Path, default=None,
        help="Curated taxonomy overrides TSV (default: DATA_DIR/plant_taxonomy_overrides.tsv).",
    )
    parser.add_argument(
        "--offline", action="store_true",
        help="Do not query GBIF; unresolved names remain Unresolved.",
    )
    parser.add_argument("--taxonomy-workers", type=int, default=8)
    parser.add_argument(
        "--refresh-taxonomy", action="store_true",
        help="Re-query GBIF for all plant names and replace cached matches.",
    )
    parser.add_argument(
        "--fasta", action="append", type=Path, default=[],
        help="FASTA file for SH extraction; repeatable. Defaults to FASTA files in data-dir.",
    )
    parser.add_argument(
        "--only", action="append", default=[],
        help="Process only input filenames containing this text; repeatable.",
    )
    parser.add_argument(
        "--prepare-taxonomy-only", action="store_true",
        help="Resolve and save taxonomy, then stop.",
    )
    args = parser.parse_args()
    if min(
        args.fungi_rank_max, args.fungi_min_count,
        args.plant_rank_max, args.plant_min_count,
    ) < 1:
        parser.error("All rank and count thresholds must be positive integers")
    if args.n_randomizations < 2:
        parser.error("--n-randomizations must be at least 2")
    if args.jobs < 0:
        parser.error("--jobs must be 0 (automatic) or a positive integer")
    if not math.isfinite(args.figure_font_scale) or args.figure_font_scale <= 0:
        parser.error("--figure-font-scale must be a finite number greater than 0")
    args.analysis_level = analysis_level
    return args


def is_unknown(value: object) -> bool:
    return str(value).strip().casefold() in UNKNOWN_LABELS


def split_sample_id_tokens(value: object) -> list[str]:
    """Split sample_id/sample_ids fields into stable, de-duplicated tokens."""
    if value is None:
        return []
    if isinstance(value, float) and np.isnan(value):
        return []
    text = str(value).strip()
    if not text or text.casefold() in UNKNOWN_LABELS:
        return []
    tokens = [
        token.strip()
        for token in re.split(r"[;,|]+", text)
        if token.strip() and token.strip().casefold() not in UNKNOWN_LABELS
    ]
    return sorted(set(tokens))


def combine_sample_id_values(values: Iterable[object]) -> str:
    """Combine possibly delimited sample identifiers with a delimiter R can parse."""
    tokens: set[str] = set()
    for value in values:
        tokens.update(split_sample_id_tokens(value))
    return "|".join(sorted(tokens))


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("_") or "matrix"


def read_count_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig", index_col=0)
    if df.empty:
        raise ValueError(f"Empty input matrix: {path}")
    if df.index.has_duplicates:
        duplicates = df.index[df.index.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Duplicate fungal row names in {path}: {duplicates}")
    if df.columns.duplicated().any():
        duplicates = df.columns[df.columns.duplicated()].tolist()[:5]
        raise ValueError(f"Duplicate plant columns in {path}: {duplicates}")
    numeric = df.apply(pd.to_numeric, errors="coerce")
    bad = numeric.isna() & ~df.isna()
    if bad.any().any():
        row, col = np.argwhere(bad.to_numpy())[0]
        raise ValueError(
            f"Non-numeric count in {path}, fungal row {df.index[row]!r}, plant column {df.columns[col]!r}."
        )
    numeric = numeric.fillna(0.0)
    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"Counts must be finite and non-negative: {path}")
    if not np.allclose(values, np.rint(values)):
        raise ValueError(f"Null models require integer counts: {path}")
    numeric = numeric.astype(np.int64)
    numeric.index = numeric.index.astype(str).str.strip()
    numeric.columns = numeric.columns.astype(str).str.strip()
    return numeric


def remove_unknown_and_empty(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    row_keep = np.array([not is_unknown(x) for x in df.index])
    col_keep = np.array([not is_unknown(x) for x in df.columns])
    out = df.loc[row_keep, col_keep].copy()
    empty_rows = out.sum(axis=1).eq(0)
    empty_cols = out.sum(axis=0).eq(0)
    out = out.loc[~empty_rows, ~empty_cols]
    report = {
        "unknown_rows_removed": int((~row_keep).sum()),
        "unknown_columns_removed": int((~col_keep).sum()),
        "zero_rows_removed_after_unknown_filter": int(empty_rows.sum()),
        "zero_columns_removed_after_unknown_filter": int(empty_cols.sum()),
        "remaining_rows": int(out.shape[0]),
        "remaining_columns": int(out.shape[1]),
        "remaining_total_count": int(out.to_numpy().sum()),
    }
    if out.empty or out.to_numpy().sum() == 0:
        raise ValueError("No non-zero data remain after removing Unknown and empty rows/columns.")
    return out, report


def sort_matrix(df: pd.DataFrame) -> pd.DataFrame:
    row_totals = df.sum(axis=1)
    row_order = sorted(df.index, key=lambda x: (-int(row_totals[x]), str(x).casefold(), str(x)))
    col_totals = df.sum(axis=0)
    col_order = sorted(df.columns, key=lambda x: (-int(col_totals[x]), str(x).casefold(), str(x)))
    return df.loc[row_order, col_order]


def write_tsv(df: pd.DataFrame, path: Path, index: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    clean = df.copy()
    for column in clean.select_dtypes(include=["object", "string"]).columns:
        clean[column] = clean[column].map(
            lambda x: re.sub(r"[\t\r\n\u2028\u2029]+", " ", x) if isinstance(x, str) else x
        )
    clean.to_csv(path, sep="\t", index=index, encoding="utf-8", lineterminator="\n", na_rep="NA")
    verify_tsv(path, expected_rows=len(clean), expected_columns=len(clean.columns) + int(index))


def verify_tsv(path: Path, expected_rows: int, expected_columns: int) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) != expected_columns:
            raise RuntimeError(f"Invalid TSV header in {path}")
        count = 0
        for number, row in enumerate(reader, start=2):
            count += 1
            if len(row) != expected_columns:
                raise RuntimeError(f"TSV column mismatch in {path}, record {number}")
            if any(len(cell) > EXCEL_CELL_LIMIT for cell in row):
                raise RuntimeError(f"Excel cell limit exceeded in {path}, record {number}")
        if count != expected_rows:
            raise RuntimeError(f"TSV row mismatch in {path}: expected {expected_rows}, found {count}")


TAXONOMY_COLUMNS = [
    "input_name", "accepted_name", "accepted_rank", "genus", "family", "order",
    "class", "phylum", "kingdom", "taxon_key", "match_type", "confidence", "issues",
    "source_url", "resolved_at_utc",
]


def empty_taxonomy_record(name: str, source_url: str = "") -> dict[str, object]:
    record: dict[str, object] = {column: "" for column in TAXONOMY_COLUMNS}
    record.update({"input_name": name, "source_url": source_url, "match_type": "NONE"})
    return record


def gbif_request(name: str, force_genus: bool) -> tuple[dict[str, object], str]:
    params = {"name": name, "kingdom": "Plantae"}
    if force_genus:
        params["rank"] = "GENUS"
    url = GBIF_MATCH_URL + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "specificity-network-workflow/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.loads(response.read().decode("utf-8")), url


def resolve_one_plant_name(name: str) -> dict[str, object]:
    last_url = ""
    for force_genus in (True, False):
        try:
            payload, last_url = gbif_request(name, force_genus)
        except Exception as exc:
            record = empty_taxonomy_record(name, last_url)
            record["issues"] = f"request_error:{type(exc).__name__}:{exc}"
            return record
        match_type = str(payload.get("matchType", "NONE")).upper()
        matched_rank = str(payload.get("rank", "")).upper()
        # A rank-constrained query can return only Plantae/Tracheophyta as a
        # HIGHERRANK match. Retry without the GENUS constraint so that family,
        # order, tribe-like, and legacy names can receive a more useful match.
        if force_genus and matched_rank != "GENUS":
            continue
        if match_type != "NONE":
            record = empty_taxonomy_record(name, last_url)
            mapping = {
                "scientificName": "accepted_name", "rank": "accepted_rank", "genus": "genus",
                "family": "family", "order": "order", "class": "class", "phylum": "phylum",
                "kingdom": "kingdom", "usageKey": "taxon_key", "matchType": "match_type",
                "confidence": "confidence",
            }
            for source, target in mapping.items():
                value = payload.get(source, "")
                record[target] = "" if value is None else value
            issues = payload.get("issues", [])
            record["issues"] = ";".join(issues) if isinstance(issues, list) else str(issues or "")
            record["resolved_at_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
            return record
    return empty_taxonomy_record(name, last_url)


def load_taxonomy(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=TAXONOMY_COLUMNS)
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig", dtype=str, keep_default_na=False)
    if "input_name" not in df.columns:
        raise ValueError(f"Taxonomy file needs an input_name column: {path}")
    for column in TAXONOMY_COLUMNS:
        if column not in df:
            df[column] = ""
    return df[TAXONOMY_COLUMNS].drop_duplicates("input_name", keep="last")


def prepare_taxonomy(
    plant_names: Iterable[str], path: Path, offline: bool, workers: int,
    refresh: bool = False, override_path: Path | None = None,
) -> pd.DataFrame:
    names = sorted(set(str(x).strip() for x in plant_names if not is_unknown(x)), key=str.casefold)
    cached = load_taxonomy(path)
    original_cached = cached.copy()
    cached_names = set() if refresh else set(cached["input_name"])
    missing = [name for name in names if name not in cached_names]
    if missing and offline:
        additions = pd.DataFrame([empty_taxonomy_record(name) for name in missing])
    elif missing:
        print(f"Resolving {len(missing)} plant name(s) with the GBIF species-match API...")
        records: list[dict[str, object]] = []
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            futures = {pool.submit(resolve_one_plant_name, name): name for name in missing}
            for completed, future in enumerate(as_completed(futures), start=1):
                records.append(future.result())
                if completed % 50 == 0 or completed == len(futures):
                    print(f"  taxonomy: {completed}/{len(futures)}")
        additions = pd.DataFrame(records)
    else:
        additions = pd.DataFrame(columns=TAXONOMY_COLUMNS)
    if refresh:
        cached = cached.loc[~cached["input_name"].isin(names)]
    combined = pd.concat([cached, additions], ignore_index=True)
    combined = combined.drop_duplicates("input_name", keep="last")
    if override_path is not None and override_path.exists():
        overrides = pd.read_csv(
            override_path, sep="\t", encoding="utf-8-sig", dtype=str, keep_default_na=False
        )
        if "input_name" not in overrides:
            raise ValueError(f"Taxonomy override file needs input_name: {override_path}")
        combined = combined.set_index("input_name", drop=False)
        for _, override in overrides.iterrows():
            name = override["input_name"]
            if name not in combined.index:
                combined.loc[name] = empty_taxonomy_record(name)
            for column in TAXONOMY_COLUMNS:
                if column in override and str(override[column]).strip():
                    combined.loc[name, column] = override[column]
        combined = combined.reset_index(drop=True)
    combined = combined[TAXONOMY_COLUMNS].sort_values("input_name", key=lambda s: s.str.casefold())
    original_sorted = original_cached[TAXONOMY_COLUMNS].sort_values(
        "input_name", key=lambda s: s.str.casefold()
    ).reset_index(drop=True)
    combined_reset = combined.reset_index(drop=True)
    # Avoid concurrent rewrites of an unchanged shared taxonomy cache when
    # several independent input groups are analysed in parallel.
    if not original_sorted.equals(combined_reset):
        write_tsv(combined_reset, path, index=False)
    return combined_reset


def taxonomy_for_names(taxonomy: pd.DataFrame, names: Iterable[str]) -> pd.DataFrame:
    wanted = list(dict.fromkeys(str(x) for x in names))
    indexed = taxonomy.set_index("input_name", drop=False)
    records = []
    for name in wanted:
        if name in indexed.index:
            records.append(indexed.loc[name].to_dict())
        else:
            records.append(empty_taxonomy_record(name))
    return pd.DataFrame(records, columns=TAXONOMY_COLUMNS)


def remove_invalid_plant_columns(
    analysis: pd.DataFrame, taxonomy: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Remove columns that are not resolved as plant genera before analysis."""
    tax = taxonomy_for_names(taxonomy, analysis.columns).set_index("input_name", drop=False)
    invalid_genus_labels = {
        "viridiplantae", "plantae", "tracheophyta", "streptophyta",
        "embryophyta", "vascular plants", "plant", "plants",
    }
    keep: list[str] = []
    excluded: list[dict[str, object]] = []
    for name in analysis.columns:
        record = tax.loc[name] if name in tax.index else pd.Series(dtype=object)
        accepted_rank = str(record.get("accepted_rank", "")).strip().upper()
        genus = str(record.get("genus", "")).strip()
        reasons: list[str] = []
        if accepted_rank != "GENUS":
            reasons.append(f"accepted_rank={accepted_rank or 'unresolved'}")
        if not genus:
            reasons.append("genus_missing")
        elif genus.casefold() in invalid_genus_labels:
            reasons.append(f"invalid_genus={genus}")
        if reasons:
            excluded.append({
                "input_plant_name": name,
                "removed_total_count": int(analysis[name].sum()),
                "removal_reason": ";".join(reasons),
                "accepted_name": record.get("accepted_name", ""),
                "accepted_rank": accepted_rank,
                "genus": genus,
                "family": record.get("family", ""),
                "order": record.get("order", ""),
                "source_url": record.get("source_url", ""),
            })
        else:
            keep.append(name)
    out = analysis.loc[:, keep].copy()
    zero_rows = out.sum(axis=1).eq(0)
    out = out.loc[~zero_rows]
    if out.empty or out.to_numpy().sum() == 0:
        raise ValueError("No non-zero data remain after removing invalid plant-genus columns.")
    excluded_df = pd.DataFrame(excluded, columns=[
        "input_plant_name", "removed_total_count", "removal_reason", "accepted_name",
        "accepted_rank", "genus", "family", "order", "source_url",
    ])
    report = {
        "invalid_plant_columns_removed": len(excluded),
        "counts_removed_with_invalid_plant_columns": int(
            sum(item["removed_total_count"] for item in excluded)
        ),
        "zero_fungal_rows_removed_after_invalid_plant_filter": int(zero_rows.sum()),
    }
    return out, excluded_df, report


def rank_label(record: pd.Series, input_name: str, rank: str) -> str:
    if rank == "genus":
        value = str(record.get("genus", "")).strip()
        return value or input_name
    value = str(record.get(rank, "")).strip()
    accepted_rank = str(record.get("accepted_rank", "")).strip().casefold()
    accepted_name = str(record.get("accepted_name", "")).strip()
    if not value and accepted_rank == rank:
        value = accepted_name
    return value or "Unresolved"


def aggregate_plant_rank(
    analysis: pd.DataFrame, taxonomy: pd.DataFrame, rank: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tax = taxonomy_for_names(taxonomy, analysis.columns)
    tax_index = tax.set_index("input_name")
    mapping = {
        name: rank_label(tax_index.loc[name], name, rank) if name in tax_index.index else "Unresolved"
        for name in analysis.columns
    }
    grouped = analysis.T.groupby(pd.Index([mapping[x] for x in analysis.columns]), sort=False).sum().T
    grouped.index.name = analysis.index.name or "fungus"
    grouped.columns.name = f"plant_{rank}"
    grouped = sort_matrix(grouped)
    tax.insert(1, "aggregation_rank", rank)
    tax.insert(2, "aggregation_label", [mapping[x] for x in tax["input_name"]])
    return grouped, tax


def select_matrix(
    df: pd.DataFrame,
    fungi_selection: str,
    fungi_rank_max: int,
    fungi_min_count: int,
    plant_selection: str,
    plant_rank_max: int,
    plant_min_count: int,
    fungal_sampling_unit_counts: pd.Series,
    plant_sampling_unit_counts: pd.Series,
) -> pd.DataFrame:
    """Select both axes by sampling-unit presence rank or count threshold."""
    row_totals = fungal_sampling_unit_counts.reindex(df.index, fill_value=0).astype(np.int64)
    row_order = sorted(
        df.index, key=lambda x: (-int(row_totals[x]), str(x).casefold(), str(x))
    )
    row_names = (
        list(row_order[:fungi_rank_max])
        if fungi_selection == "rank"
        else [label for label in row_order if row_totals[label] >= fungi_min_count]
    )
    eligible_columns = [x for x in df.columns if not is_unknown(x) and x != "Unresolved"]
    column_totals = plant_sampling_unit_counts.reindex(
        eligible_columns, fill_value=0
    ).astype(np.int64)
    column_order = sorted(
        eligible_columns,
        key=lambda x: (-int(column_totals[x]), str(x).casefold(), str(x)),
    )
    col_names = (
        column_order[:plant_rank_max]
        if plant_selection == "rank"
        else [label for label in column_order if column_totals[label] >= plant_min_count]
    )
    if not row_names:
        raise ValueError("No fungi passed the requested rank/count threshold.")
    if not col_names:
        raise ValueError("No plant lineages passed the requested rank/count threshold.")
    selected = df.loc[row_names, col_names]
    selected = selected.loc[selected.sum(axis=1).gt(0), selected.sum(axis=0).gt(0)]
    if selected.empty:
        raise ValueError("The selected matrix has no non-zero interactions.")
    return selected


def sampling_unit_rank_counts(
    occurrences: pd.DataFrame,
    taxonomy_for_file: pd.DataFrame,
    matrix: pd.DataFrame,
) -> tuple[pd.Series, pd.Series]:
    """Count distinct reconstructed sampling units for every fungal and plant label."""
    fungal_counts = (
        occurrences.drop_duplicates(["sampling_unit_id", "fungus"])
        .groupby("fungus")["sampling_unit_id"].nunique()
        .reindex(matrix.index, fill_value=0).astype(np.int64)
    )
    lineage = dict(zip(
        taxonomy_for_file["input_name"].astype(str),
        taxonomy_for_file["aggregation_label"].astype(str),
    ))
    plant_units = occurrences[["sampling_unit_id", "plant_genus"]].drop_duplicates()
    plant_units["plant_lineage"] = plant_units["plant_genus"].map(lineage).fillna("Unresolved")
    plant_counts = (
        plant_units.groupby("plant_lineage")["sampling_unit_id"].nunique()
        .reindex(matrix.columns, fill_value=0).astype(np.int64)
    )
    return fungal_counts, plant_counts


def sampling_unit_ranking_table(
    counts: pd.Series,
    label_column: str,
    selected_labels: Iterable[str],
    selection_mode: str,
    rank_max: int,
    min_count: int,
    ineligible_labels: Iterable[str] = (),
) -> pd.DataFrame:
    """Create an auditable sampling-unit ranking and selection-label table."""
    ineligible = {str(value) for value in ineligible_labels}
    order = sorted(
        counts.index,
        key=lambda x: (-int(counts[x]), str(x).casefold(), str(x)),
    )
    selected = {str(value) for value in selected_labels}
    records = []
    eligible_rank = 0
    for rank, label in enumerate(order, start=1):
        eligible = str(label) not in ineligible and not is_unknown(label)
        if eligible:
            eligible_rank += 1
        passes = eligible and (
            eligible_rank <= rank_max
            if selection_mode == "rank"
            else int(counts[label]) >= min_count
        )
        records.append({
            "rank": rank,
            "eligible_rank": eligible_rank if eligible else pd.NA,
            label_column: label,
            "n_sampling_units": int(counts[label]),
            "eligible_for_randomization": eligible,
            "passes_requested_threshold": passes,
            "randomization_selection": (
                "selected" if str(label) in selected else "not_selected"
            ),
        })
    return pd.DataFrame(records)


def plot_count_heatmap(df: pd.DataFrame, path: Path, title: str) -> None:
    values = df.to_numpy(dtype=float)
    masked = np.ma.masked_where(values <= 0, values)
    maximum = float(values.max())
    width = max(9.0, 0.30 * df.shape[1] + 4)
    height = max(8.0, 0.20 * df.shape[0] + 3)
    fig, ax = plt.subplots(figsize=(width, height))
    if maximum > 1:
        image = ax.imshow(masked, aspect="auto", cmap="viridis", norm=LogNorm(vmin=1, vmax=maximum))
        label = "Interaction count (log color scale)"
    else:
        image = ax.imshow(masked, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        label = "Interaction count"
    image.cmap.set_bad("white")
    ax.set_xticks(np.arange(df.shape[1]), labels=df.columns, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(np.arange(df.shape[0]), labels=df.index, fontsize=7)
    ax.set_xlabel(df.columns.name or "Host plant")
    ax.set_ylabel(df.index.name or "Fungus")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, pad=0.01, label=label)
    apply_figure_font_scale(fig)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def plot_bipartite_network(df: pd.DataFrame, path: Path, title: str) -> None:
    values = df.to_numpy(dtype=float)
    row_totals = values.sum(axis=1)
    col_totals = values.sum(axis=0)
    row_y = np.linspace(1, 0, len(df.index))
    col_y = np.linspace(1, 0, len(df.columns))
    fig_height = max(10.0, 0.22 * max(len(df.index), len(df.columns)) + 3)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    max_count = max(1.0, values.max())
    for i, j in zip(*np.nonzero(values)):
        count = values[i, j]
        width = 0.25 + 3.0 * math.sqrt(count / max_count)
        ax.plot([0, 1], [row_y[i], col_y[j]], color="#667788", alpha=0.15, lw=width, zorder=1)
    row_sizes = 25 + 220 * np.sqrt(row_totals / max(1.0, row_totals.max()))
    col_sizes = 25 + 220 * np.sqrt(col_totals / max(1.0, col_totals.max()))
    ax.scatter(np.zeros_like(row_y), row_y, s=row_sizes, color="#7B3294", edgecolor="white", zorder=3)
    ax.scatter(np.ones_like(col_y), col_y, s=col_sizes, color="#008837", edgecolor="white", zorder=3)
    for y, label in zip(row_y, df.index):
        ax.text(-0.025, y, str(label), ha="right", va="center", fontsize=7)
    for y, label in zip(col_y, df.columns):
        ax.text(1.025, y, str(label), ha="left", va="center", fontsize=8)
    ax.text(0, 1.035, "Fungi", ha="center", va="bottom", weight="bold")
    ax.text(1, 1.035, df.columns.name or "Host plants", ha="center", va="bottom", weight="bold")
    ax.set_xlim(-0.42, 1.42)
    ax.set_ylim(-0.03, 1.08)
    ax.set_title(title)
    ax.axis("off")
    apply_figure_font_scale(fig)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def z_score(observed: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    result = np.full(np.shape(observed), np.nan, dtype=float)
    valid = np.isfinite(observed) & np.isfinite(mean) & np.isfinite(sd) & (sd > 0)
    result[valid] = (observed[valid] - mean[valid]) / sd[valid]
    return result


def normal_two_sided_p(z: np.ndarray) -> np.ndarray:
    values = np.asarray(z, dtype=float)
    out = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    out[finite] = np.array([math.erfc(abs(x) / math.sqrt(2.0)) for x in values[finite]])
    return out


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan)
    flat = p.ravel()
    finite_indices = np.flatnonzero(np.isfinite(flat))
    if finite_indices.size == 0:
        return out
    order = finite_indices[np.argsort(flat[finite_indices])]
    ranked = flat[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out.ravel()[order] = np.clip(adjusted, 0, 1)
    return out


@dataclass
class NullResult:
    model: str
    observed_row_dprime: np.ndarray
    row_dprime_mean: np.ndarray
    row_dprime_sd: np.ndarray
    row_dprime_n: np.ndarray
    observed_col_dprime: np.ndarray
    col_dprime_mean: np.ndarray
    col_dprime_sd: np.ndarray
    col_dprime_n: np.ndarray
    cell_mean: np.ndarray
    cell_sd: np.ndarray
    cell_n: np.ndarray
    row_p_randomization_upper: np.ndarray
    row_p_randomization_lower: np.ndarray
    row_p_randomization_directional_tail: np.ndarray
    row_p_randomization_two_sided: np.ndarray
    col_p_randomization_upper: np.ndarray
    col_p_randomization_lower: np.ndarray
    col_p_randomization_directional_tail: np.ndarray
    col_p_randomization_two_sided: np.ndarray
    cell_p_randomization_two_sided: np.ndarray

    @property
    def row_z(self) -> np.ndarray:
        return z_score(self.observed_row_dprime, self.row_dprime_mean, self.row_dprime_sd)

    @property
    def col_z(self) -> np.ndarray:
        return z_score(self.observed_col_dprime, self.col_dprime_mean, self.col_dprime_sd)

    def cell_z(self, observed: np.ndarray) -> np.ndarray:
        return z_score(observed, self.cell_mean, self.cell_sd)


def stable_seed(base_seed: int, *parts: str) -> int:
    digest = hashlib.sha256("\0".join(parts).encode("utf-8")).digest()
    return (base_seed + int.from_bytes(digest[:8], "little")) % (2**63 - 1)


def dprime_table(
    labels: Iterable[str], observed: np.ndarray, mean: np.ndarray, sd: np.ndarray,
    n_valid: np.ndarray, p_upper: np.ndarray, p_lower: np.ndarray,
    p_directional: np.ndarray, p_two_sided: np.ndarray, label_column: str,
) -> pd.DataFrame:
    z = z_score(observed, mean, sd)
    p = normal_two_sided_p(z)
    return pd.DataFrame({
        label_column: list(labels),
        "dprime_observed": observed,
        "dprime_null_mean": mean,
        "dprime_null_sd": sd,
        "z_standardized_dprime": z,
        "p_randomization_upper": p_upper,
        "p_randomization_lower": p_lower,
        "p_randomization_directional_tail": p_directional,
        "p_randomization_two_sided": p_two_sided,
        "fdr_bh_randomization_two_sided": bh_fdr(p_two_sided),
        "p_normal_two_sided": p,
        "fdr_bh": bh_fdr(p),
        "n_valid_randomizations": n_valid,
    })


def matrix_from_array(array: np.ndarray, template: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(array, index=template.index, columns=template.columns)


def plot_preference_figure(
    selected: pd.DataFrame,
    result: NullResult,
    path: Path,
    title: str,
    column_role: str = "plant",
) -> None:
    cell_z = result.cell_z(selected.to_numpy(dtype=float))
    finite = np.abs(cell_z[np.isfinite(cell_z)])
    limit = max(3.0, float(np.nanpercentile(finite, 98)) if finite.size else 3.0)
    limit = min(limit, 10.0)
    actual_max = float(finite.max()) if finite.size else float("nan")
    fig = plt.figure(figsize=(max(11, 0.32 * selected.shape[1] + 6), max(9, 0.21 * selected.shape[0] + 4)))
    grid = fig.add_gridspec(
        2, 2,
        width_ratios=[1.15, 8],
        height_ratios=[0.7, 8],
        wspace=0.08,
        hspace=0.05,
    )
    top = fig.add_subplot(grid[0, 1])
    left = fig.add_subplot(grid[1, 0])
    heat = fig.add_subplot(grid[1, 1])
    x = np.arange(selected.shape[1])
    y = np.arange(selected.shape[0])
    top.bar(x, result.col_z, color=np.where(result.col_z >= 0, "#B2182B", "#2166AC"), width=0.85)
    top.axhline(0, color="black", lw=0.7)
    top.axhline(1.96, color="gray", lw=0.6, ls="--")
    top.axhline(-1.96, color="gray", lw=0.6, ls="--")
    top.set_xlim(-0.5, len(x) - 0.5)
    top.set_ylabel(f"{column_role.capitalize()}\nz-d'")
    top.tick_params(axis="x", labelbottom=False)
    left.barh(y, result.row_z, color=np.where(result.row_z >= 0, "#B2182B", "#2166AC"), height=0.85)
    left.axvline(0, color="black", lw=0.7)
    left.axvline(1.96, color="gray", lw=0.6, ls="--")
    left.axvline(-1.96, color="gray", lw=0.6, ls="--")
    left.set_ylim(len(y) - 0.5, -0.5)
    left.set_xlabel("Fungal z-d'")
    left.set_yticks(y, labels=selected.index, fontsize=7)
    image = heat.imshow(
        cell_z, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    )
    heat.set_xticks(x, labels=selected.columns, rotation=55, ha="right", fontsize=8)
    heat.set_yticks(y)
    heat.tick_params(axis="y", labelleft=False)
    heat.set_xlabel(selected.columns.name or column_role.capitalize())
    heat.set_ylabel("")
    colorbar = fig.colorbar(image, ax=[top, heat], fraction=0.025, pad=0.02, extend="both")
    colorbar.set_label(
        f"z-standardized 2DP (display ±{limit:.2g}; max |z|={actual_max:.2g})"
    )
    fig.suptitle(title, y=0.995)
    save_figure(fig, path)
    plt.close(fig)


def plot_z_fdr_relationships(
    selected: pd.DataFrame,
    result: NullResult,
    path: Path,
    title: str,
    column_role: str = "plant",
) -> None:
    """Draw Fig. 4B-style relationships between z scores and BH-FDR."""
    cell_z = result.cell_z(selected.to_numpy(dtype=float))
    row_fdr = bh_fdr(result.row_p_randomization_two_sided)
    col_fdr = bh_fdr(result.col_p_randomization_two_sided)
    cell_fdr = bh_fdr(result.cell_p_randomization_two_sided)
    panels = [
        (result.row_z, row_fdr, "Fungal d'", "#7B3294"),
        (result.col_z, col_fdr, f"{column_role.capitalize()} d'", "#008837"),
        (cell_z.ravel(), cell_fdr.ravel(), "Two-dimensional preference (2DP)", "#444444"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    for ax, (z_values, fdr_values, label, color) in zip(axes, panels):
        z_values = np.asarray(z_values, dtype=float).ravel()
        fdr_values = np.asarray(fdr_values, dtype=float).ravel()
        finite = np.isfinite(z_values) & np.isfinite(fdr_values)
        ax.scatter(
            z_values[finite], fdr_values[finite], s=16, alpha=0.55,
            color=color, edgecolors="none",
        )
        ax.axhline(
            0.05, color="#B2182B", lw=1.0, ls="--",
            label="two-sided FDR = 0.05 (0.025 per tail)",
        )
        ax.axvline(0, color="gray", lw=0.7)
        ax.axvline(-3, color="gray", lw=0.7, ls=":")
        ax.axvline(3, color="gray", lw=0.7, ls=":")
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel(f"z-standardized {label}")
        ax.set_title(label)
        ax.grid(alpha=0.15)
    axes[0].set_ylabel("Benjamini-Hochberg FDR")
    axes[-1].legend(loc="upper right", frameon=False)
    fig.suptitle(title)
    apply_figure_font_scale(fig)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def write_randomization_outputs(
    selected: pd.DataFrame,
    result: NullResult,
    output_dir: Path,
    prefix: str,
    n: int,
    column_role: str = "plant",
) -> None:
    model = result.model
    row_table = dprime_table(
        selected.index, result.observed_row_dprime, result.row_dprime_mean,
        result.row_dprime_sd, result.row_dprime_n,
        result.row_p_randomization_upper, result.row_p_randomization_lower,
        result.row_p_randomization_directional_tail,
        result.row_p_randomization_two_sided,
        selected.index.name or "fungus",
    )
    col_table = dprime_table(
        selected.columns, result.observed_col_dprime, result.col_dprime_mean,
        result.col_dprime_sd, result.col_dprime_n,
        result.col_p_randomization_upper, result.col_p_randomization_lower,
        result.col_p_randomization_directional_tail,
        result.col_p_randomization_two_sided,
        selected.columns.name or "plant_lineage",
    )
    write_tsv(row_table, output_dir / f"{prefix}_{model}_fungal_dprime.tsv", index=False)
    write_tsv(
        col_table,
        output_dir / f"{prefix}_{model}_{safe_name(column_role)}_dprime.tsv",
        index=False,
    )
    observed = selected.to_numpy(dtype=float)
    cell_z = result.cell_z(observed)
    cell_p = normal_two_sided_p(cell_z)
    cell_fdr_randomization = bh_fdr(result.cell_p_randomization_two_sided)
    for label, array in [
        ("2dp_z", cell_z), ("2dp_p_normal", cell_p), ("2dp_fdr_bh", bh_fdr(cell_p)),
        ("2dp_p_randomization_two_sided", result.cell_p_randomization_two_sided),
        ("2dp_fdr_bh_randomization_two_sided", cell_fdr_randomization),
        ("null_cell_mean", result.cell_mean), ("null_cell_sd", result.cell_sd),
    ]:
        write_tsv(matrix_from_array(array, selected), output_dir / f"{prefix}_{model}_{label}.tsv")
    row_index, col_index = np.indices(selected.shape)
    cell_long = pd.DataFrame({
        selected.index.name or "fungus": np.asarray(selected.index)[row_index.ravel()],
        selected.columns.name or "plant_lineage": np.asarray(selected.columns)[col_index.ravel()],
        "observed_count": observed.ravel(),
        "z_standardized_2dp": cell_z.ravel(),
        "p_randomization_two_sided": result.cell_p_randomization_two_sided.ravel(),
        "fdr_bh_randomization_two_sided": cell_fdr_randomization.ravel(),
    })
    write_tsv(
        cell_long,
        output_dir / f"{prefix}_{model}_2dp_z_fdr_long.tsv",
        index=False,
    )
    plot_preference_figure(
        selected, result, output_dir / f"{prefix}_{model}_z_dprime_2dp.pdf",
        f"z-standardized d' and 2DP ({model}; {n:,} randomizations)",
        column_role,
    )
    plot_z_fdr_relationships(
        selected, result, output_dir / f"{prefix}_{model}_z_fdr_relationships.pdf",
        f"z-standardized preference and FDR ({model}; {n:,} randomizations)",
        column_role,
    )


def fasta_paths(data_dir: Path, explicit: list[Path]) -> list[Path]:
    if explicit:
        return explicit
    paths: list[Path] = []
    for pattern in ("*.fasta", "*.fa", "*.fna", "*.fasta.gz", "*.fa.gz", "*.fna.gz"):
        paths.extend(data_dir.glob(pattern))
    return sorted(set(paths))


def open_text(path: Path):
    if path.suffix == ".gz":
        import gzip
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def iter_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    identifier = re.split(r"[|\s]", header[1:], maxsplit=1)[0]
                    yield identifier, header, "".join(chunks)
                header = line
                chunks = []
            elif header is not None:
                chunks.append(line)
        if header is not None:
            identifier = re.split(r"[|\s]", header[1:], maxsplit=1)[0]
            yield identifier, header, "".join(chunks)


def extract_selected_fasta(
    selected_ids: list[str], paths: list[Path], output_fasta: Path, missing_tsv: Path,
    confidence_label: str,
) -> None:
    wanted = set(selected_ids)
    found: dict[str, tuple[str, str]] = {}
    for path in paths:
        for identifier, header, sequence in iter_fasta(path):
            if identifier in wanted and identifier not in found:
                found[identifier] = (header, sequence)
        if len(found) == len(wanted):
            break
    with output_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier in selected_ids:
            if identifier not in found:
                continue
            header, sequence = found[identifier]
            handle.write(header + "\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    report = pd.DataFrame({
        "sh_id": selected_ids,
        "fasta_found": [identifier in found for identifier in selected_ids],
        "confidence_subset": confidence_label,
        "fasta_source_files": ";".join(path.name for path in paths),
    })
    write_tsv(report, missing_tsv, index=False)
    print(f"    FASTA: {len(found)}/{len(selected_ids)} selected SH sequences found")


def is_sh_matrix(df: pd.DataFrame) -> bool:
    return bool(len(df.index)) and all(re.fullmatch(r"SH\d+(?:\.\d+FU)?", x) for x in df.index)


def confidence_label_from_name(name: str) -> Optional[str]:
    """In v2, only the medium-or-higher subset is analyzed."""
    normalized = name.casefold().replace("-", "_")
    if re.search(r"confidence_high(?:\.|_)", normalized):
        return None
    if re.search(r"confidence_medium_or_high(?:\.|_)", normalized):
        return "medium_or_high"
    if "root_occurrences_with_hosts" in normalized or "root_host_sequences" in normalized:
        return "medium_or_high"
    return None


def fasta_paths_for_confidence(
    paths: list[Path], confidence_label: str
) -> list[Path]:
    """Select only FASTA files whose header counts match the input confidence."""
    if not paths:
        return []
    matched = [
        path for path in paths
        if confidence_label_from_name(path.name) == confidence_label
    ]
    if not matched:
        raise FileNotFoundError(
            f"No FASTA file for confidence_{confidence_label} was found. "
            "Rerun 1_GlobalFungi to create confidence-specific FASTA files or "
            "supply them with --fasta. Available FASTA files: "
            + ", ".join(path.name for path in paths)
        )
    def its_priority(path: Path) -> tuple[int, str]:
        name = path.name.casefold()
        if "its1_and_its2" in name:
            return (0, name)
        if re.search(r"_all(?:\.|_)", name):
            return (1, name)
        if "its1_only" in name:
            return (2, name)
        if "its2_only" in name:
            return (3, name)
        return (4, name)

    return sorted(matched, key=its_priority)


def input_paths(
    data_dir: Path, pattern: str, excluded_paths: Iterable[Path], only: list[str]
) -> list[Path]:
    excluded = {path.resolve() for path in excluded_paths}
    paths = []
    for path in sorted(data_dir.glob(pattern)):
        if path.resolve() in excluded:
            continue
        if only and not any(token.casefold() in path.name.casefold() for token in only):
            continue
        confidence = confidence_label_from_name(path.name)
        if confidence not in CONFIDENCE_LABELS:
            continue
        paths.append(path)
    if not paths:
        raise FileNotFoundError(f"No input TSV files matched {data_dir / pattern}")
    return paths


def threshold_token(selection: str, rank_max: int, min_count: int) -> str:
    return str(rank_max) if selection == "rank" else f"min{min_count}"


def analysis_run_name(args: argparse.Namespace) -> str:
    fungi_label = "sh" if args.analysis_level == "sh" else "genus"
    fungi_token = f"{fungi_label}{threshold_token(args.fungi_selection, args.fungi_rank_max, args.fungi_min_count)}"
    plant_token = f"{args.plant_rank}{threshold_token(args.plant_selection, args.plant_rank_max, args.plant_min_count)}"
    return f"{fungi_token}_{plant_token}_{args.n_randomizations}"


def resolve_user_path(path: Path, project_dir: Path) -> Path:
    candidates = [path]
    if not path.is_absolute():
        candidates.extend([Path.cwd() / path, project_dir / path, project_dir.parent / path])
    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved.exists():
            return resolved
    return path.expanduser().resolve()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable_report_path(path: Union[Path, str], *bases: Path) -> str:
    candidate = Path(path)
    try:
        resolved = candidate.resolve()
    except OSError:
        return str(candidate)
    for base in bases:
        try:
            return str(resolved.relative_to(base.resolve()))
        except ValueError:
            continue
    return f"[EXTERNAL]/{resolved.name}"



def v2_medium_or_high_filename(source_name: str) -> str:
    """Use neutral filenames inside v2 runs while sourcing medium-or-higher inputs."""
    name = source_name
    name = name.replace("_confidence_medium_or_high", "")
    name = name.replace("confidence_medium_or_high_", "")
    name = name.replace("confidence_medium_or_high", "")
    name = re.sub(r"__+", "_", name)
    return name

def copy_globalfungi_inputs(args: argparse.Namespace, project_dir: Path) -> None:
    if args.globalfungi_results is None:
        return
    source_dir = resolve_user_path(args.globalfungi_results, project_dir)
    if not source_dir.is_dir():
        raise FileNotFoundError(f"--globalfungi-results is not a directory: {source_dir}")
    run_name = analysis_run_name(args)
    output_root = args.output_dir
    if output_root.name != run_name:
        output_root = output_root / run_name
    args.output_dir = output_root.resolve()
    copied_data_dir = args.output_dir / "data"
    copied_data_dir.mkdir(parents=True, exist_ok=True)

    required_pattern_groups = [
        ["*root_occurrences_with_hosts.tsv", "*root_occurrences_with_hosts_confidence_medium_or_high.tsv"],
    ]
    optional_patterns = [
        "*root_host_sequences_all.fasta",
        "*root_host_sequences_its1_and_its2.fasta",
        "*root_host_sequences_its1_only.fasta",
        "*root_host_sequences_its2_only.fasta",
        "*root_host_sequence_summary_all.tsv",
        "*root_host_sequence_summary_its1_and_its2.tsv",
        "*root_host_sequence_summary_its1_only.tsv",
        "*root_host_sequence_summary_its2_only.tsv",
        "*root_host_sequences_confidence_medium_or_high*.fasta",
        "*root_host_sequence_summary_confidence_medium_or_high*.tsv",
        "*fasta_linkage_report.tsv",
        "*unique_occurrence_count_validation.tsv",
        "dataset_summary.tsv",
        "run_configuration.tsv",
        "run_command.sh",
        "run_timing.tsv",
        "runtime_versions.tsv",
    ]
    rows: list[dict[str, object]] = []
    copied: set[Path] = set()
    for patterns in required_pattern_groups:
        matches: list[Path] = []
        for pattern in patterns:
            matches = sorted(source_dir.glob(pattern))
            if matches:
                break
        if not matches:
            raise FileNotFoundError(
                "Required GlobalFungi input not found; tried: "
                + ", ".join(str(source_dir / pattern) for pattern in patterns)
            )
        for source in matches:
            copied.add(source.resolve())
            destination = copied_data_dir / v2_medium_or_high_filename(source.name)
            shutil.copy2(source, destination)
            rows.append({
                "role": "required_occurrence",
                "source_file": portable_report_path(source, project_dir.parent),
                "copied_file": str(destination.relative_to(args.output_dir)),
                "bytes": destination.stat().st_size,
                "sha256": file_sha256(destination),
            })
    for pattern in optional_patterns:
        for source in sorted(source_dir.glob(pattern)):
            if source.resolve() in copied:
                continue
            copied.add(source.resolve())
            destination = copied_data_dir / v2_medium_or_high_filename(source.name)
            shutil.copy2(source, destination)
            rows.append({
                "role": "supporting_input",
                "source_file": portable_report_path(source, project_dir.parent),
                "copied_file": str(destination.relative_to(args.output_dir)),
                "bytes": destination.stat().st_size,
                "sha256": file_sha256(destination),
            })
    if args.taxonomy is None:
        for source in (project_dir / "data" / "plant_taxonomy.tsv",
                       project_dir / "data" / "plant_taxonomy_overrides.tsv"):
            if source.exists():
                destination = copied_data_dir / v2_medium_or_high_filename(source.name)
                shutil.copy2(source, destination)
                rows.append({
                    "role": "local_taxonomy_cache",
                    "source_file": portable_report_path(source, project_dir.parent),
                    "copied_file": str(destination.relative_to(args.output_dir)),
                    "bytes": destination.stat().st_size,
                    "sha256": file_sha256(destination),
                })
    write_tsv(
        pd.DataFrame(rows),
        args.output_dir / "globalfungi_input_copy_manifest.tsv",
        index=False,
    )
    args.data_dir = copied_data_dir
    args.globalfungi_results = source_dir
    args.globalfungi_results_report = portable_report_path(source_dir, project_dir.parent)
    args.analysis_run_name = run_name


def occurrence_path_for_input(path: Path, args: argparse.Namespace) -> Path:
    confidence = confidence_label_from_name(path.name)
    if confidence not in CONFIDENCE_LABELS:
        raise ValueError(f"Cannot determine the v2 medium-or-higher input subset from {path.name}")
    return args.occurrence_paths[confidence]


def runtime_version_records(args: argparse.Namespace) -> list[dict[str, str]]:
    """Record the full Python/R environments plus scripts and run configuration."""
    records: list[dict[str, str]] = []

    def add(category: str, name: object, version: object = "", path: object = "", details: object = "") -> None:
        records.append({
            "category": str(category), "name": str(name), "version": str(version),
            "path": str(path), "details": str(details),
        })

    add("run", "timestamp", datetime.now(timezone.utc).astimezone().isoformat())
    add("run", "command", details=" ".join([sys.executable, *sys.argv]))
    add("run", "analysis_level", args.analysis_level)
    add("system", "platform", platform.platform())
    add("system", "machine", platform.machine())
    add("python", platform.python_implementation(), platform.python_version(), sys.executable, sys.version)
    script_paths = {
        Path(sys.argv[0]).resolve(), Path(__file__).resolve(), args.r_script.resolve()
    }
    for script_path in sorted(script_paths):
        if script_path.exists():
            add("script", script_path.name, hashlib.sha256(script_path.read_bytes()).hexdigest(), script_path, "SHA-256")
    distributions = []
    for distribution in importlib_metadata.distributions():
        distributions.append((
            distribution.metadata.get("Name") or "Unknown",
            distribution.version,
            distribution.locate_file(""),
        ))
    for name, version, location in sorted(distributions, key=lambda x: (str(x[0]).casefold(), str(x[1]))):
        add("python_package", name, version, location)
    rscript = shutil.which("Rscript")
    if rscript is None:
        add("R", "Rscript", "not_found", details="Rscript is required for randomization analysis")
        return records
    version = subprocess.run([rscript, "--version"], capture_output=True, text=True, check=False)
    add("R", "R", (version.stdout or version.stderr).strip(), rscript)
    packages = subprocess.run(
        [
            rscript, "-e",
            "ip <- installed.packages()[,c('Package','Version','LibPath'),drop=FALSE]; "
            "write.table(ip,sep='\\t',row.names=FALSE,col.names=FALSE,quote=FALSE)",
        ],
        capture_output=True, text=True, check=False,
    )
    if packages.returncode == 0:
        for line in packages.stdout.splitlines():
            fields = line.split("\t", 2)
            if len(fields) >= 2:
                add("R_package", fields[0], fields[1], fields[2] if len(fields) == 3 else "")
    else:
        add("R", "installed_packages", "inspection_failed", details=packages.stderr.strip())
    return records


def append_log(path: Path, message: str) -> None:
    timestamp = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"[{timestamp}] {message}\n")


def write_run_configuration(args: argparse.Namespace, inputs: list[Path]) -> None:
    """Write a reusable shell command, resolved settings, and an initial run log."""
    command = [sys.executable, *sys.argv]
    command_text = shlex.join(command)
    command_path = args.output_dir / "run_command.sh"
    with command_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        handle.write(command_text + "\n")
    configuration = {
        "command": command_text,
        "analysis_level": args.analysis_level,
        "input_files": ";".join(portable_report_path(path, args.output_dir) for path in inputs),
        "analysis_run_name": getattr(args, "analysis_run_name", ""),
        "globalfungi_results": getattr(args, "globalfungi_results_report", ""),
        "copied_input_data_dir": portable_report_path(args.data_dir, args.output_dir),
        "input_subset": "medium_or_higher_only",
        "occurrence_input_glob": args.input_glob,
        "plant_rank": args.plant_rank,
        "fungi_selection": args.fungi_selection,
        "fungal_ranking_metric": "number of sampling units with fungal presence",
        "fungi_rank_max": args.fungi_rank_max,
        "fungi_min_count": args.fungi_min_count,
        "plant_selection": args.plant_selection,
        "plant_ranking_metric": "number of sampling units assigned to plant lineage",
        "plant_rank_max": args.plant_rank_max,
        "plant_min_count": args.plant_min_count,
        "n_randomizations": args.n_randomizations,
        "figure_font_scale": args.figure_font_scale,
        "randomization_method": (
            "global_sampling_unit_plant_label_shuffle;"
            "global_sampling_unit_continent_label_shuffle"
        ),
        "seed": args.seed,
        "jobs_requested": args.jobs,
        "input_workers": args.effective_input_jobs,
        "total_worker_budget": args.effective_jobs,
        "occurrence_unit": "binary fungal presence per latitude-longitude-host sampling unit",
        "randomization_unit": "latitude-longitude-host sampling unit",
        "randomization_constraint": (
            "shuffle sampling-unit plant or continent labels globally without replacement; "
            "no blocks"
        ),
        "matrix_reconstruction": (
            "rebuild fungus-by-plant or fungus-by-continent matrix after every shuffle"
        ),
        "continent_specificity_analysis": True,
        "two_sided_tests": True,
        "fdr_method": "Benjamini-Hochberg on empirical two-sided P values",
    }
    write_tsv(
        pd.DataFrame(
            [{"setting": key, "value": value} for key, value in configuration.items()]
        ),
        args.output_dir / "run_configuration.tsv",
        index=False,
    )
    log_path = args.output_dir / "workflow.log"
    log_path.write_text("", encoding="utf-8")
    append_log(log_path, f"Started: {command_text}")
    append_log(log_path, f"Timer start: {RUN_STATE.get('started_at', 'unknown')}")
    append_log(log_path, f"Inputs: {configuration['input_files']}")
    append_log(
        log_path,
        "Occurrence unit: binary fungal presence per latitude-longitude-host sampling unit",
    )
    RUN_STATE.update({
        "output_dir": args.output_dir,
        "log_path": log_path,
        "command": command_text,
    })


def finalize_run_timing(status: str, exit_code: int) -> None:
    """Append end time and elapsed wall-clock time, including failed runs."""
    if "log_path" not in RUN_STATE or RUN_STATE.get("finalized"):
        return
    ended_at = datetime.now(timezone.utc).astimezone()
    elapsed = time.perf_counter() - float(RUN_STATE["started_monotonic"])
    append_log(Path(RUN_STATE["log_path"]), f"End: status={status}, exit_code={exit_code}")
    append_log(Path(RUN_STATE["log_path"]), f"Elapsed seconds: {elapsed:.6f}")
    write_tsv(
        pd.DataFrame([{
            "start_time": RUN_STATE["started_at"],
            "end_time": ended_at.isoformat(),
            "elapsed_seconds": elapsed,
            "status": status,
            "exit_code": exit_code,
            "command": RUN_STATE["command"],
        }]),
        Path(RUN_STATE["output_dir"]) / "run_timing.tsv",
        index=False,
    )
    RUN_STATE["finalized"] = True


def require_r_environment(args: argparse.Namespace) -> None:
    if shutil.which("Rscript") is None:
        raise RuntimeError("Rscript was not found on PATH.")
    if not args.r_script.exists():
        raise FileNotFoundError(f"R occurrence-shuffle helper does not exist: {args.r_script}")
    required = ["bipartite"]
    expression = (
        "missing <- c(" + ",".join(repr(name) for name in required) + ")"
        "[!vapply(c(" + ",".join(repr(name) for name in required) + "), requireNamespace, logical(1), quietly=TRUE)];"
        "if(length(missing)){cat(paste(missing,collapse=','));quit(status=2)}"
    )
    check = subprocess.run(["Rscript", "-e", expression], capture_output=True, text=True, check=False)
    if check.returncode != 0:
        missing = check.stdout.strip() or check.stderr.strip()
        raise RuntimeError(
            f"Required R package(s) are unavailable: {missing}. Run "
            "Rscript -e \"install.packages('bipartite')\"."
        )


def _cell_array(cell_table: pd.DataFrame, selected: pd.DataFrame, value: str) -> np.ndarray:
    matrix = cell_table.pivot(index="row_label", columns="col_label", values=value)
    return matrix.reindex(index=selected.index, columns=selected.columns).to_numpy(dtype=float)


def run_r_sample_unit_shuffle(
    selected: pd.DataFrame,
    occurrence_labels: pd.DataFrame,
    n_randomizations: int,
    seed: int,
    r_script: Path,
    output_dir: Path,
    prefix: str,
    write_outputs: bool = True,
    column_role: str = "plant",
) -> tuple[NullResult, pd.DataFrame]:
    """Shuffle one sampling-unit label globally without spatial blocks."""
    output_dir.mkdir(parents=True, exist_ok=True)
    model = "global_sampling_unit_label_shuffle"
    token = safe_name(f"{prefix}_{model}")
    input_path = output_dir / f".{token}__r_input.tsv"
    r_prefix = output_dir / f".{token}__r_result"
    occurrence_labels.to_csv(
        input_path, sep="\t", index=False, encoding="utf-8", lineterminator="\n"
    )
    r_seed = int(seed % 2_147_483_646) + 1
    command = [
        "Rscript", str(r_script), str(input_path), str(r_prefix),
        str(n_randomizations), str(r_seed), column_role,
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    r_log_path = output_dir / f"{prefix}_{model}_R.log"
    with r_log_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("COMMAND\n" + shlex.join(command) + "\n\nSTDOUT\n")
        handle.write(completed.stdout)
        handle.write("\nSTDERR\n")
        handle.write(completed.stderr)
    if completed.stdout.strip():
        print(completed.stdout.strip())
    if completed.stderr.strip():
        print(completed.stderr.strip())
    if completed.returncode != 0:
        input_path.unlink(missing_ok=True)
        for suffix in ["__row.tsv", "__col.tsv", "__cell.tsv", "__constraints.tsv"]:
            Path(f"{r_prefix}{suffix}").unlink(missing_ok=True)
        raise RuntimeError(f"R occurrence-label shuffle failed: {completed.stderr.strip()}")
    row_table = pd.read_csv(f"{r_prefix}__row.tsv", sep="\t", encoding="utf-8")
    col_table = pd.read_csv(f"{r_prefix}__col.tsv", sep="\t", encoding="utf-8")
    cell_table = pd.read_csv(f"{r_prefix}__cell.tsv", sep="\t", encoding="utf-8")
    constraints = pd.read_csv(f"{r_prefix}__constraints.tsv", sep="\t", encoding="utf-8")
    if row_table["label"].astype(str).tolist() != selected.index.astype(str).tolist():
        raise RuntimeError("R row labels changed during occurrence-label shuffling.")
    if col_table["label"].astype(str).tolist() != selected.columns.astype(str).tolist():
        raise RuntimeError("R column labels changed during occurrence-label shuffling.")
    result = NullResult(
        model=model,
        observed_row_dprime=row_table["dprime_observed"].to_numpy(float),
        row_dprime_mean=row_table["dprime_null_mean"].to_numpy(float),
        row_dprime_sd=row_table["dprime_null_sd"].to_numpy(float),
        row_dprime_n=row_table["n_valid_randomizations"].to_numpy(np.int64),
        observed_col_dprime=col_table["dprime_observed"].to_numpy(float),
        col_dprime_mean=col_table["dprime_null_mean"].to_numpy(float),
        col_dprime_sd=col_table["dprime_null_sd"].to_numpy(float),
        col_dprime_n=col_table["n_valid_randomizations"].to_numpy(np.int64),
        cell_mean=_cell_array(cell_table, selected, "null_mean"),
        cell_sd=_cell_array(cell_table, selected, "null_sd"),
        cell_n=_cell_array(cell_table, selected, "n_valid_randomizations").astype(np.int64),
        row_p_randomization_upper=row_table["p_randomization_upper"].to_numpy(float),
        row_p_randomization_lower=row_table["p_randomization_lower"].to_numpy(float),
        row_p_randomization_directional_tail=row_table[
            "p_randomization_directional_tail"
        ].to_numpy(float),
        row_p_randomization_two_sided=row_table[
            "p_randomization_two_sided"
        ].to_numpy(float),
        col_p_randomization_upper=col_table["p_randomization_upper"].to_numpy(float),
        col_p_randomization_lower=col_table["p_randomization_lower"].to_numpy(float),
        col_p_randomization_directional_tail=col_table[
            "p_randomization_directional_tail"
        ].to_numpy(float),
        col_p_randomization_two_sided=col_table[
            "p_randomization_two_sided"
        ].to_numpy(float),
        cell_p_randomization_two_sided=_cell_array(cell_table, selected, "p_randomization_two_sided"),
    )
    row_label = selected.index.name or "fungus"
    row_table = row_table.rename(columns={"label": row_label})
    write_tsv(
        constraints,
        output_dir / f"{prefix}_{model}_randomization_constraints.tsv",
        index=False,
    )
    if write_outputs:
        write_randomization_outputs(
            selected,
            result,
            output_dir,
            prefix,
            n_randomizations,
            column_role,
        )
    for suffix in ["__row.tsv", "__col.tsv", "__cell.tsv", "__constraints.tsv"]:
        Path(f"{r_prefix}{suffix}").unlink(missing_ok=True)
    input_path.unlink(missing_ok=True)
    return result, row_table


OCCURRENCE_ALIASES = {
    "sample_id": ("sample_id", "sample_ID"),
    "sh_id": ("sh_id", "sequence_id"),
    "genus": ("genus", "fungal_genus"),
    "continent": ("continent",),
    "country": ("country",),
    "latitude": ("latitude",),
    "longitude": ("longitude",),
    "ph": ("ph", "pH"),
    "mat": ("mat", "MAT_study"),
    "map": ("map", "MAP_study"),
    "soil_carbon": ("soil_carbon", "organic_C_content"),
    "species": ("species",), "family": ("family",), "order": ("order",),
    "class": ("class",), "phylum": ("phylum",), "kingdom": ("kingdom",),
    "taxonomy": ("taxonomy",),
    "source_study": ("source_study", "study_title", "doi"),
    "host_candidate": ("host_candidate",),
    "occurrence_count": ("occurrence_count",),
    "unique_occurrence_count": ("unique_occurrence_count",),
    "n_sample_sh_occurrences": ("n_sample_sh_occurrences",),
    "sample_ids": ("sample_ids",),
    "source_studies": ("source_studies",),
}


def occurrence_column_map(path: Path) -> dict[str, str]:
    columns = pd.read_csv(path, sep="\t", nrows=0, encoding="utf-8-sig").columns.tolist()
    lookup = {column.casefold(): column for column in columns}
    mapping = {}
    for standard, aliases in OCCURRENCE_ALIASES.items():
        for alias in aliases:
            if alias.casefold() in lookup:
                mapping[standard] = lookup[alias.casefold()]
                break
    required = {
        "sample_id", "sh_id", "continent", "latitude", "longitude",
        "host_candidate", "unique_occurrence_count",
    }
    missing = required - mapping.keys()
    if missing:
        raise ValueError(f"Occurrence table is missing required columns: {sorted(missing)}")
    return mapping


def host_plant_genus(value: object) -> str:
    """Extract one plausible plant genus; ambiguous/higher-rank values are Unknown."""
    if pd.isna(value):
        return "Unknown"
    text = str(value).strip()
    if is_unknown(text) or re.search(r"[,;/|]", text) or re.search(
        r"\b(?:and|mixed|multiple|various|several)\b", text, flags=re.I
    ):
        return "Unknown"
    first = text.split()[0].strip("()[]{}.,;:'\"")
    if not re.fullmatch(r"[A-Z][A-Za-z-]+", first):
        return "Unknown"
    if first.casefold().endswith(("aceae", "oideae", "phyta", "opsida")):
        return "Unknown"
    return first


def read_spatial_occurrences(
    path: Path, analysis_level: str, block_column: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read and validate unique occurrences, returning kept and excluded rows."""
    mapping = occurrence_column_map(path)
    header = pd.read_csv(path, sep="\t", nrows=0, encoding="utf-8-sig").columns
    lookup = {column.casefold(): column for column in header}
    if block_column.casefold() not in lookup:
        raise ValueError(f"Spatial-block column {block_column!r} is absent from {path}")
    usecols = sorted(set(mapping.values()) | {lookup[block_column.casefold()]})
    data = pd.read_csv(path, sep="\t", usecols=usecols, encoding="utf-8-sig", low_memory=False)
    data = data.rename(columns={source: standard for standard, source in mapping.items()})
    block_source = lookup[block_column.casefold()]
    if block_source in data.columns and block_source != "spatial_block":
        data["spatial_block"] = data[block_source]
    elif block_column in data.columns:
        data["spatial_block"] = data[block_column]
    fungus_column = "sh_id" if analysis_level == "sh" else "genus"
    data["fungus"] = data[fungus_column].astype("string").str.strip()
    data["plant_genus"] = data["host_candidate"].map(host_plant_genus)
    data["spatial_block"] = data["spatial_block"].astype("string").str.strip()
    data["unique_occurrence_count"] = pd.to_numeric(
        data["unique_occurrence_count"], errors="coerce"
    )
    reasons = pd.Series("", index=data.index, dtype="string")

    def add_reason(mask: pd.Series, reason: str) -> None:
        current = reasons.loc[mask].fillna("")
        reasons.loc[mask] = np.where(current.eq(""), reason, current + ";" + reason)

    add_reason(data["unique_occurrence_count"].ne(1) | data["unique_occurrence_count"].isna(), "invalid_unique_occurrence_count")
    add_reason(data["fungus"].isna() | data["fungus"].map(is_unknown), "unknown_fungus")
    if analysis_level == "sh":
        add_reason(~data["fungus"].str.fullmatch(r"SH\d+(?:\.\d+FU)?", na=False), "invalid_sh_id")
    add_reason(data["plant_genus"].map(is_unknown), "unknown_or_ambiguous_host_genus")
    add_reason(data["spatial_block"].isna() | data["spatial_block"].map(is_unknown), "unknown_spatial_block")
    latitude = pd.to_numeric(data["latitude"], errors="coerce")
    longitude = pd.to_numeric(data["longitude"], errors="coerce")
    add_reason(latitude.isna() | longitude.isna(), "missing_coordinates")
    complete = (
        latitude.notna()
        & longitude.notna()
        & data["host_candidate"].notna()
    )
    duplicated = complete & data.duplicated(
        ["sh_id", "latitude", "longitude", "host_candidate"], keep=False
    )
    add_reason(duplicated, "duplicated_unique_occurrence_key")
    excluded = data.loc[reasons.ne("")].copy()
    excluded.insert(0, "exclusion_reason", reasons.loc[excluded.index])
    kept = data.loc[reasons.eq("")].copy()
    kept["occurrence_count"] = 1
    kept["latitude"] = latitude.loc[kept.index]
    kept["longitude"] = longitude.loc[kept.index]
    sampling_keys = pd.MultiIndex.from_frame(
        kept[["spatial_block", "latitude", "longitude", "host_candidate"]]
    )
    kept["sampling_unit_id"] = pd.factorize(sampling_keys, sort=False)[0] + 1
    return kept.reset_index(drop=True), excluded.reset_index(drop=True)


def occurrence_matrix_by_genus(occurrences: pd.DataFrame) -> pd.DataFrame:
    binary = occurrences.drop_duplicates(
        ["sampling_unit_id", "fungus", "plant_genus"]
    )
    matrix = pd.crosstab(
        binary["fungus"], binary["plant_genus"],
        values=binary["occurrence_count"], aggfunc="sum",
    ).fillna(0).astype(np.int64)
    matrix.index.name = "fungus"
    matrix.columns.name = "plant_genus"
    return sort_matrix(matrix)


def sampling_unit_continent_matrix(
    occurrences: pd.DataFrame, selected_fungi: Iterable[str]
) -> pd.DataFrame:
    """Count sampling-unit fungal presences by continent for specificity analysis."""
    binary = occurrences.drop_duplicates(
        ["sampling_unit_id", "fungus", "spatial_block"]
    )
    matrix = pd.crosstab(binary["fungus"], binary["spatial_block"])
    matrix = matrix.reindex(index=list(selected_fungi), fill_value=0).astype(np.int64)
    matrix = matrix.loc[matrix.sum(axis=1).gt(0), matrix.sum(axis=0).gt(0)]
    matrix.index.name = occurrences.attrs.get("fungus_label", "fungus")
    matrix.columns.name = "continent"
    column_order = matrix.sum(axis=0).sort_values(ascending=False).index
    return matrix.loc[:, column_order]


def sample_level_binary_matrix(
    occurrences: pd.DataFrame, fungi_order: Iterable[str]
) -> pd.DataFrame:
    """Return sampling units x fungi with one binary presence per cell."""
    binary = occurrences.drop_duplicates(["sampling_unit_id", "fungus"])
    matrix = pd.crosstab(binary["sampling_unit_id"], binary["fungus"])
    matrix = matrix.reindex(columns=list(fungi_order), fill_value=0).astype(np.int8)
    matrix = matrix.loc[matrix.sum(axis=1).gt(0)]
    matrix.index.name = "sampling_unit_id"
    if not np.isin(matrix.to_numpy(), [0, 1]).all():
        raise RuntimeError("Sample-level matrix contains a value other than 0 or 1.")
    return matrix


def sampling_unit_metadata(
    occurrences: pd.DataFrame,
    taxonomy_for_file: pd.DataFrame,
) -> pd.DataFrame:
    """Describe each reconstructed site-host sampling unit."""
    lineage = dict(zip(
        taxonomy_for_file["input_name"].astype(str),
        taxonomy_for_file["aggregation_label"].astype(str),
    ))
    work = occurrences.copy()
    work["plant_lineage"] = work["plant_genus"].map(lineage).fillna("Unresolved")
    invariant = [
        "spatial_block", "continent", "latitude", "longitude",
        "host_candidate", "plant_genus", "plant_lineage", "country",
    ]
    invariant = [column for column in invariant if column in work.columns]
    for column in invariant:
        inconsistent = work.groupby("sampling_unit_id", sort=False)[column].nunique(
            dropna=True
        ).gt(1)
        if inconsistent.any():
            examples = inconsistent[inconsistent].index.astype(str).tolist()[:5]
            raise RuntimeError(
                f"Sampling-unit metadata column {column!r} is inconsistent for: {examples}"
            )
    metadata = (
        work.sort_values("sampling_unit_id")
        .drop_duplicates("sampling_unit_id")
        .set_index("sampling_unit_id")[invariant]
    )

    def joined_values(series: pd.Series) -> str:
        output: list[str] = []
        seen: set[str] = set()
        for value in series.dropna().astype(str):
            for item in (part.strip() for part in value.split(";")):
                if item and not is_unknown(item) and item not in seen:
                    seen.add(item)
                    output.append(item)
        return "; ".join(output)

    for column in ["sample_id", "sample_ids", "source_study", "source_studies"]:
        if column in work.columns:
            metadata[column] = work.groupby("sampling_unit_id", sort=False)[column].agg(
                joined_values
            ).reindex(metadata.index)
    metadata["n_binary_fungal_presences"] = (
        work.drop_duplicates(["sampling_unit_id", "fungus"])
        .groupby("sampling_unit_id").size().reindex(metadata.index).astype(np.int64)
    )
    return metadata.reset_index()


def selected_occurrence_label_table(
    occurrences: pd.DataFrame,
    taxonomy_for_file: pd.DataFrame,
    selected: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = dict(zip(
        taxonomy_for_file["input_name"].astype(str),
        taxonomy_for_file["aggregation_label"].astype(str),
    ))
    work = occurrences.copy()
    work["plant_lineage"] = work["plant_genus"].map(mapping).fillna("Unresolved")
    work = work[
        work["fungus"].isin(selected.index.astype(str))
        & work["plant_lineage"].isin(selected.columns.astype(str))
    ]
    sample_columns = [column for column in ("sample_id", "sample_ids") if column in work.columns]
    if sample_columns:
        work["_shuffle_sample_ids"] = work[sample_columns].apply(
            lambda row: combine_sample_id_values(row.values), axis=1
        )
    else:
        work["_shuffle_sample_ids"] = ""
    occurrence_labels = (
        work.groupby(
            ["sampling_unit_id", "spatial_block", "fungus", "plant_lineage"],
            sort=False,
            dropna=False,
        )
        .agg(sample_ids=("_shuffle_sample_ids", combine_sample_id_values))
        .reset_index()
    )
    occurrence_labels = occurrence_labels.rename(columns={
        "fungus": "row_label", "plant_lineage": "col_label",
    })
    row_order = {label: i for i, label in enumerate(selected.index.astype(str), start=1)}
    col_order = {label: i for i, label in enumerate(selected.columns.astype(str), start=1)}
    occurrence_labels["row_order"] = occurrence_labels["row_label"].map(row_order)
    occurrence_labels["col_order"] = occurrence_labels["col_label"].map(col_order)
    occurrence_labels.insert(
        0, "occurrence_id", np.arange(1, len(occurrence_labels) + 1, dtype=np.int64)
    )
    block_summary = (
        work.groupby("spatial_block", sort=False)
        .agg(
            n_binary_fungus_occurrences=("occurrence_count", "sum"),
            n_sampling_units=("sampling_unit_id", "nunique"),
            n_fungi=("fungus", "nunique"),
            n_plant_lineages=("plant_lineage", "nunique"),
        )
        .reset_index()
    )
    if len(occurrence_labels) != int(selected.to_numpy().sum()):
        raise RuntimeError("Occurrence rows do not equal the selected host matrix total.")
    reconstructed = pd.crosstab(
        occurrence_labels["row_label"], occurrence_labels["col_label"]
    ).reindex(index=selected.index, columns=selected.columns, fill_value=0)
    if not np.array_equal(
        reconstructed.to_numpy(dtype=np.int64), selected.to_numpy(dtype=np.int64)
    ):
        raise RuntimeError("Occurrence labels do not reconstruct the selected host matrix.")
    return occurrence_labels, block_summary


def continent_occurrence_label_table(
    occurrences: pd.DataFrame,
    selected: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the binary sampling-unit table used to shuffle continent labels."""
    work = occurrences[
        occurrences["fungus"].isin(selected.index.astype(str))
        & occurrences["spatial_block"].isin(selected.columns.astype(str))
    ].copy()
    work = work.drop_duplicates(["sampling_unit_id", "fungus", "spatial_block"])
    occurrence_labels = work[[
        "sampling_unit_id", "spatial_block", "fungus"
    ]].copy()
    occurrence_labels["col_label"] = occurrence_labels["spatial_block"]
    occurrence_labels = occurrence_labels.rename(columns={"fungus": "row_label"})
    row_order = {label: i for i, label in enumerate(selected.index.astype(str), start=1)}
    col_order = {label: i for i, label in enumerate(selected.columns.astype(str), start=1)}
    occurrence_labels["row_order"] = occurrence_labels["row_label"].map(row_order)
    occurrence_labels["col_order"] = occurrence_labels["col_label"].map(col_order)
    occurrence_labels.insert(
        0, "occurrence_id", np.arange(1, len(occurrence_labels) + 1, dtype=np.int64)
    )
    summary = (
        work.groupby("spatial_block", sort=False)
        .agg(
            n_binary_fungus_occurrences=("occurrence_count", "sum"),
            n_sampling_units=("sampling_unit_id", "nunique"),
            n_fungi=("fungus", "nunique"),
        )
        .reset_index()
        .rename(columns={"spatial_block": "continent"})
    )
    if len(occurrence_labels) != int(selected.to_numpy().sum()):
        raise RuntimeError("Occurrence rows do not equal the selected continent matrix total.")
    reconstructed = pd.crosstab(
        occurrence_labels["row_label"], occurrence_labels["col_label"]
    ).reindex(index=selected.index, columns=selected.columns, fill_value=0)
    if not np.array_equal(
        reconstructed.to_numpy(dtype=np.int64), selected.to_numpy(dtype=np.int64)
    ):
        raise RuntimeError("Occurrence labels do not reconstruct the continent matrix.")
    return occurrence_labels, summary


def load_selected_occurrences(
    path: Path, selected_fungi: Iterable[str], analysis_level: str
) -> pd.DataFrame:
    """Stream unique occurrences without re-collapsing them by sample ID."""
    mapping = occurrence_column_map(path)
    fungus_column = "sh_id" if analysis_level == "sh" else "genus"
    if fungus_column not in mapping:
        raise ValueError(f"Occurrence table lacks the {fungus_column} column.")
    selected_set = {str(value) for value in selected_fungi}
    usecols = sorted(set(mapping.values()))
    chunks = []
    for chunk in pd.read_csv(
        path, sep="\t", usecols=usecols, chunksize=100_000,
        encoding="utf-8-sig", low_memory=False,
    ):
        rename = {source: standard for standard, source in mapping.items()}
        chunk = chunk.rename(columns=rename)
        chunk["fungus"] = chunk[fungus_column].astype("string").str.strip()
        chunk["fungus"] = chunk["fungus"].fillna("Unknown")
        chunk.loc[chunk["fungus"].map(is_unknown), "fungus"] = "Unknown"
        chunk = chunk[chunk["fungus"].isin(selected_set)]
        if not chunk.empty:
            chunks.append(chunk)
    if not chunks:
        return pd.DataFrame(columns=["fungus", "sample_id", *mapping.keys()])
    data = pd.concat(chunks, ignore_index=True)
    for column in ["latitude", "longitude", "ph", "mat", "map", "soil_carbon"]:
        if column in data:
            data[column] = pd.to_numeric(data[column], errors="coerce")
    data["sample_id"] = data["sample_id"].astype("string").str.strip()
    if "sh_id" in data:
        data["sh_id"] = data["sh_id"].astype("string").str.strip()
    data["host_candidate"] = data["host_candidate"].astype("string").str.strip()
    data["unique_occurrence_count"] = pd.to_numeric(
        data["unique_occurrence_count"], errors="coerce"
    )
    invalid_count = data["unique_occurrence_count"].ne(1) | data[
        "unique_occurrence_count"
    ].isna()
    if invalid_count.any():
        raise ValueError(
            f"Unique-occurrence input contains {int(invalid_count.sum()):,} row(s) "
            "whose unique_occurrence_count is not 1."
        )
    complete = (
        data["latitude"].notna()
        & data["longitude"].notna()
        & data["host_candidate"].notna()
        & ~data["host_candidate"].map(is_unknown)
    )
    unique_keys = ["sh_id", "latitude", "longitude", "host_candidate"]
    duplicated = complete & data.duplicated(unique_keys, keep=False)
    if duplicated.any():
        raise ValueError(
            f"Unique-occurrence input contains {int(duplicated.sum()):,} duplicated "
            "complete SH-latitude-longitude-host keys. Regenerate it with 1_GlobalFungi."
        )
    data["occurrence_count"] = data["unique_occurrence_count"].astype(np.int64)
    data["unique_occurrence_row_id"] = np.arange(1, len(data) + 1)
    return data


def first_nonmissing(series: pd.Series) -> object:
    values = series.dropna().astype(str).str.strip()
    values = values[~values.str.casefold().isin(UNKNOWN_LABELS)]
    return values.iloc[0] if not values.empty else pd.NA


def integrated_fungus_summary(
    selected: pd.DataFrame,
    occurrences: pd.DataFrame,
    host_rows: dict[str, pd.DataFrame],
    continent_rows: dict[str, pd.DataFrame],
    analysis_level: str,
) -> pd.DataFrame:
    labels = list(selected.index.astype(str))
    base = pd.DataFrame({"fungus": labels})
    totals = selected.sum(axis=1)
    base["selection_rank"] = np.arange(1, len(base) + 1)
    base["host_matrix_total_count"] = base["fungus"].map(totals)
    taxonomic = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "taxonomy"]
    aggregation: dict[str, tuple[str, object]] = {
        "n_occurrences": ("occurrence_count", "sum"),
        "n_sample_ids": ("sample_id", "nunique"),
        "n_countries": ("country", "nunique"),
        "n_continents": ("continent", "nunique"),
    }
    for variable in ["latitude", "ph", "mat", "map", "soil_carbon"]:
        if variable in occurrences:
            aggregation[f"median_{variable}"] = (variable, "median")
    for column in taxonomic:
        if column in occurrences:
            aggregation[column] = (column, first_nonmissing)
    if occurrences.empty:
        environmental = pd.DataFrame({"fungus": labels})
    else:
        environmental = occurrences.groupby("fungus", sort=False).agg(**aggregation).reset_index()
    summary = base.merge(environmental, on="fungus", how="left")
    for domain, tables in (("host", host_rows), ("continent", continent_rows)):
        for model, table in tables.items():
            label_column = table.columns[0]
            subset = table[[
                label_column, "dprime_observed", "dprime_null_mean", "dprime_null_sd",
                "z_standardized_dprime", "p_randomization_two_sided",
                "fdr_bh_randomization_two_sided",
            ]].copy()
            subset = subset.rename(columns={
                label_column: "fungus",
                **{
                    column: f"{domain}_{column}_{model}"
                    for column in subset.columns[1:]
                },
            })
            summary = summary.merge(subset, on="fungus", how="left")
    output_label = "sh_id" if analysis_level == "sh" else "fungal_genus"
    return summary.rename(columns={"fungus": output_label})


def draw_world_boundaries(ax, geojson_path: Path) -> None:
    if not geojson_path.exists():
        raise FileNotFoundError(
            f"World boundary GeoJSON does not exist: {geojson_path}. "
            "Land outlines are required for spatial plots."
        )
    with geojson_path.open("r", encoding="utf-8") as handle:
        world = json.load(handle)
    for feature in world.get("features", []):
        geometry = feature.get("geometry") or {}
        coordinates = geometry.get("coordinates", [])
        polygons = [coordinates] if geometry.get("type") == "Polygon" else coordinates
        for polygon in polygons:
            for ring in polygon:
                points = np.asarray(ring, dtype=float)
                if points.ndim != 2 or points.shape[1] < 2:
                    continue
                jumps = np.flatnonzero(np.abs(np.diff(points[:, 0])) > 180) + 1
                for segment in np.split(points, jumps):
                    if len(segment) > 1:
                        ax.plot(segment[:, 0], segment[:, 1], color="#666666", lw=0.35, zorder=1)


def plot_selected_spatial_distributions(
    occurrences: pd.DataFrame,
    selected_ids: list[str],
    geojson_path: Path,
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for rank, fungus in enumerate(selected_ids, start=1):
        fungus_occurrences = occurrences.loc[occurrences["fungus"].eq(fungus)]
        n_occurrences = int(fungus_occurrences["occurrence_count"].sum())
        points = fungus_occurrences.loc[
            fungus_occurrences["latitude"].between(-90, 90)
            & fungus_occurrences["longitude"].between(-180, 180),
            ["sample_id", "latitude", "longitude", "occurrence_count"],
        ]
        fig, ax = plt.subplots(figsize=(11, 5.8))
        draw_world_boundaries(ax, geojson_path)
        if not points.empty:
            ax.scatter(
                points["longitude"], points["latitude"], s=13, alpha=0.65,
                color="#8B2F67", edgecolors="none", zorder=2,
            )
        ax.set(xlim=(-180, 180), ylim=(-90, 90), xlabel="Longitude", ylabel="Latitude")
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
        ax.grid(alpha=0.18, lw=0.4)
        ax.set_title(
            f"{fungus}: sampling-unit fungal presences "
            f"(n = {n_occurrences:,}; mapped = {int(points['occurrence_count'].sum()):,})"
        )
        apply_figure_font_scale(fig)
        fig.tight_layout()
        save_figure(
            fig,
            output_dir / f"{rank:03d}_{safe_name(fungus)}_spatial_distribution.pdf",
        )
        plt.close(fig)


def process_one(
    path: Path, taxonomy: pd.DataFrame, fasta_files: list[Path], args: argparse.Namespace
) -> None:
    # Required for multiprocessing with the spawn start method (for example macOS).
    set_figure_font_scale(args.figure_font_scale)
    input_confidence = confidence_label_from_name(path.name)
    confidence_fasta_files = fasta_files
    if args.analysis_level == "sh":
        if input_confidence is None:
            raise ValueError(f"Cannot determine confidence subset from {path.name}.")
        # Validate before the expensive randomization rather than failing at export.
        confidence_fasta_files = fasta_paths_for_confidence(
            fasta_files, input_confidence
        )
    stem = safe_name(path.stem)
    output_dir = args.output_dir / stem
    output_dir.mkdir(parents=True, exist_ok=True)
    analysis_log = output_dir / "analysis.log"
    analysis_log.write_text("", encoding="utf-8")
    append_log(analysis_log, f"Started unique-occurrence input={path}")
    print(f"\nProcessing {path.name}")
    all_occurrences, excluded_occurrences = read_spatial_occurrences(
        path, args.analysis_level, "continent"
    )
    all_occurrences.attrs["analysis_level"] = args.analysis_level
    original = occurrence_matrix_by_genus(all_occurrences)
    original.index.name = "sh_id" if args.analysis_level == "sh" else "fungal_genus"
    occurrence_row_total = int(all_occurrences["occurrence_count"].sum())
    binary_interaction_total = int(original.to_numpy().sum())
    write_tsv(
        excluded_occurrences,
        output_dir / f"{stem}_excluded_unique_occurrences.tsv",
        index=False,
    )
    write_tsv(
        pd.DataFrame([{
            "unique_occurrence_file": portable_report_path(path, args.output_dir),
            "count_unit": "binary fungal presence per latitude-longitude-host sampling unit",
            "input_rows": len(all_occurrences) + len(excluded_occurrences),
            "accepted_unique_occurrence_rows": occurrence_row_total,
            "excluded_rows": len(excluded_occurrences),
            "sampling_units": int(all_occurrences["sampling_unit_id"].nunique()),
            "binary_sampling_unit_fungus_interactions": binary_interaction_total,
            "rows_collapsed_as_duplicate_genus_presence": occurrence_row_total - binary_interaction_total,
            "matrix_total": binary_interaction_total,
            "totals_match": binary_interaction_total == int(original.to_numpy().sum()),
        }]),
        output_dir / f"{stem}_unique_occurrence_input_validation.tsv",
        index=False,
    )
    analysis, filter_report = remove_unknown_and_empty(original)
    analysis, excluded_plants, invalid_report = remove_invalid_plant_columns(analysis, taxonomy)
    filter_report.update(invalid_report)
    analysis.index.name = original.index.name or "fungus"
    analysis = sort_matrix(analysis)
    write_tsv(analysis, output_dir / f"{stem}_analysis.tsv")
    write_tsv(
        excluded_plants,
        output_dir / f"{stem}_excluded_plant_columns.tsv",
        index=False,
    )
    grouped, tax_for_file = aggregate_plant_rank(analysis, taxonomy, args.plant_rank)
    rank_prefix = f"{stem}_by_host_{args.plant_rank}"
    write_tsv(tax_for_file, output_dir / f"{stem}_plant_taxonomy.tsv", index=False)
    write_tsv(grouped, output_dir / f"{rank_prefix}.tsv")
    valid_genera = set(analysis.columns.astype(str))
    valid_fungi = set(analysis.index.astype(str))
    filtered_occurrences = all_occurrences[
        all_occurrences["plant_genus"].isin(valid_genera)
        & all_occurrences["fungus"].isin(valid_fungi)
    ].copy()
    fungal_rank_counts, plant_rank_counts = sampling_unit_rank_counts(
        filtered_occurrences, tax_for_file, grouped
    )
    selected = select_matrix(
        grouped,
        args.fungi_selection, args.fungi_rank_max, args.fungi_min_count,
        args.plant_selection, args.plant_rank_max, args.plant_min_count,
        fungal_rank_counts, plant_rank_counts,
    )
    write_tsv(selected, output_dir / f"{rank_prefix}_selected.tsv")
    fungal_label = "sh_id" if args.analysis_level == "sh" else "fungal_genus"
    fungal_ranking = sampling_unit_ranking_table(
        fungal_rank_counts, fungal_label, selected.index.astype(str),
        args.fungi_selection, args.fungi_rank_max, args.fungi_min_count,
    )
    plant_label = f"plant_{args.plant_rank}"
    plant_ranking = sampling_unit_ranking_table(
        plant_rank_counts, plant_label, selected.columns.astype(str),
        args.plant_selection, args.plant_rank_max, args.plant_min_count,
        ineligible_labels=("Unresolved", "Unknown"),
    )
    write_tsv(
        fungal_ranking,
        output_dir / f"{rank_prefix}_fungal_sampling_unit_ranking.tsv",
        index=False,
    )
    write_tsv(
        plant_ranking,
        output_dir / f"{rank_prefix}_plant_sampling_unit_ranking.tsv",
        index=False,
    )
    sample_matrix = sample_level_binary_matrix(
        filtered_occurrences, analysis.index.astype(str)
    )
    sample_metadata = sampling_unit_metadata(filtered_occurrences, tax_for_file)
    sample_metadata = sample_metadata.set_index("sampling_unit_id").reindex(
        sample_matrix.index
    ).reset_index()
    if int(sample_matrix.to_numpy().sum()) != int(analysis.to_numpy().sum()):
        raise RuntimeError(
            "Full sample-level matrix total does not equal the analysis matrix total."
        )
    write_tsv(
        sample_matrix,
        output_dir / f"{stem}_sample_level_binary_matrix.tsv",
    )
    write_tsv(
        sample_metadata,
        output_dir / f"{stem}_sample_level_metadata.tsv",
        index=False,
    )
    occurrence_labels, block_summary = selected_occurrence_label_table(
        filtered_occurrences, tax_for_file, selected
    )
    selected_presence = occurrence_labels[["sampling_unit_id", "row_label"]].rename(
        columns={"row_label": "fungus"}
    )
    selected_sample_matrix = sample_level_binary_matrix(
        selected_presence, selected.index.astype(str)
    )
    selected_sample_metadata = sample_metadata.set_index("sampling_unit_id").reindex(
        selected_sample_matrix.index
    ).reset_index()
    if int(selected_sample_matrix.to_numpy().sum()) != int(selected.to_numpy().sum()):
        raise RuntimeError(
            "Selected sample-level matrix total does not equal the selected interaction matrix total."
        )
    write_tsv(
        selected_sample_matrix,
        output_dir / f"{rank_prefix}_selected_sample_level_binary_matrix.tsv",
    )
    write_tsv(
        selected_sample_metadata,
        output_dir / f"{rank_prefix}_selected_sample_level_metadata.tsv",
        index=False,
    )
    write_tsv(
        occurrence_labels,
        output_dir / f"{rank_prefix}_selected_occurrences_for_randomization.tsv",
        index=False,
    )
    write_tsv(
        block_summary,
        output_dir / f"{rank_prefix}_selected_spatial_block_summary.tsv",
        index=False,
    )
    plot_count_heatmap(
        selected, output_dir / f"{rank_prefix}_selected_heatmap.pdf",
        f"Interaction counts: {path.stem} ({args.plant_rank} level)",
    )
    plot_bipartite_network(
        selected, output_dir / f"{rank_prefix}_selected_network.pdf",
        f"Fungus-host network: {path.stem} ({args.plant_rank} level)",
    )
    metadata = {
        "input_file": portable_report_path(path, args.output_dir),
        "analysis_level": args.analysis_level,
        "plant_rank": args.plant_rank,
        "fungi_selection": args.fungi_selection,
        "fungal_ranking_metric": "number of sampling units with fungal presence",
        "fungi_rank_max": args.fungi_rank_max,
        "fungi_min_count": args.fungi_min_count,
        "plant_selection": args.plant_selection,
        "plant_ranking_metric": "number of sampling units assigned to plant lineage",
        "plant_rank_max": args.plant_rank_max,
        "plant_min_count": args.plant_min_count,
        "selected_fungi": int(selected.shape[0]),
        "selected_plant_lineages": int(selected.shape[1]),
        "selected_total_count": int(selected.to_numpy().sum()),
        "n_randomizations": args.n_randomizations,
        "figure_font_scale": args.figure_font_scale,
        "randomization_method": (
            "global_sampling_unit_plant_label_shuffle;"
            "global_sampling_unit_continent_label_shuffle"
        ),
        "total_worker_budget": args.effective_jobs,
        "input_workers": args.effective_input_jobs,
        "seed": args.seed,
        "r_script": str(args.r_script),
        "occurrence_file": portable_report_path(path, args.output_dir),
        "count_unit": "binary fungal presence per latitude-longitude-host sampling unit",
        "randomization_constraint": (
            "sampling-unit plant and continent labels independently shuffled globally "
            "without replacement; no spatial blocks"
        ),
        "matrix_reconstruction": (
            "fungus-by-plant or fungus-by-continent matrix rebuilt after every shuffle "
            "before dprime and 2DP"
        ),
        "continent_specificity_analysis": True,
        "statistical_test": "empirical two-sided for dprime and 2DP",
        **filter_report,
    }
    write_tsv(pd.DataFrame([metadata]), output_dir / f"{stem}_analysis_metadata.tsv", index=False)

    selected_ids = list(selected.index.astype(str))
    occurrences = filtered_occurrences[
        filtered_occurrences["fungus"].isin(selected_ids)
    ].copy()
    occurrences = occurrences.drop_duplicates(["sampling_unit_id", "fungus"]).copy()
    occurrences.attrs["analysis_level"] = args.analysis_level
    occurrences.attrs["fungus_label"] = (
        "sh_id" if args.analysis_level == "sh" else "fungal_genus"
    )
    observed_continents = sampling_unit_continent_matrix(occurrences, selected_ids)
    write_tsv(
        observed_continents,
        output_dir / f"{stem}_by_continent.tsv",
    )
    model = "global_sampling_unit_label_shuffle"
    task_seed = stable_seed(args.seed, path.name, args.plant_rank, model)
    append_log(
        analysis_log,
        "Launching global sampling-unit plant-label shuffle without blocks",
    )
    print("  R global sampling-unit plant-label shuffle without blocks")
    try:
        randomization_result, row_table = run_r_sample_unit_shuffle(
            selected, occurrence_labels, args.n_randomizations, task_seed,
            args.r_script, output_dir, rank_prefix, False, "plant",
        )
    except RuntimeError as error:
        warnings.warn(str(error))
        append_log(analysis_log, f"host_global/{model}: not_estimable: {error}")
        null_status = [{
            "domain": "host_global", "model": model,
            "status": "not_estimable", "message": str(error),
        }]
        host_rows: dict[str, pd.DataFrame] = {}
    else:
        write_randomization_outputs(
            selected,
            randomization_result,
            output_dir,
            rank_prefix,
            args.n_randomizations,
            "plant",
        )
        host_rows = {model: row_table}
        null_status = [{
            "domain": "host_global", "model": model,
            "status": "completed", "message": "",
        }]
        append_log(analysis_log, f"host_global/{model}: completed")

    continent_rows: dict[str, pd.DataFrame] = {}
    continent_prefix = f"{stem}_by_continent"
    if observed_continents.shape[0] >= 2 and observed_continents.shape[1] >= 2:
        continent_labels, continent_summary = continent_occurrence_label_table(
            occurrences, observed_continents
        )
        write_tsv(
            continent_labels,
            output_dir / f"{continent_prefix}_selected_occurrences_for_randomization.tsv",
            index=False,
        )
        write_tsv(
            continent_summary,
            output_dir / f"{continent_prefix}_sampling_unit_summary.tsv",
            index=False,
        )
        continent_seed = stable_seed(args.seed, path.name, "continent", model)
        append_log(
            analysis_log,
            "Launching global sampling-unit continent-label shuffle without blocks",
        )
        print("  R global sampling-unit continent-label shuffle without blocks")
        try:
            continent_result, continent_row_table = run_r_sample_unit_shuffle(
                observed_continents,
                continent_labels,
                args.n_randomizations,
                continent_seed,
                args.r_script,
                output_dir,
                continent_prefix,
                False,
                "continent",
            )
        except RuntimeError as error:
            warnings.warn(str(error))
            append_log(analysis_log, f"continent_global/{model}: not_estimable: {error}")
            null_status.append({
                "domain": "continent_global", "model": model,
                "status": "not_estimable", "message": str(error),
            })
        else:
            write_randomization_outputs(
                observed_continents,
                continent_result,
                output_dir,
                continent_prefix,
                args.n_randomizations,
                "continent",
            )
            continent_rows = {model: continent_row_table}
            null_status.append({
                "domain": "continent_global", "model": model,
                "status": "completed", "message": "",
            })
            append_log(analysis_log, f"continent_global/{model}: completed")
    else:
        message = (
            "continent matrix requires at least two fungi and two continents; "
            f"observed dimensions={observed_continents.shape[0]}x{observed_continents.shape[1]}"
        )
        warnings.warn(message)
        append_log(analysis_log, f"continent_global/{model}: not_estimable: {message}")
        null_status.append({
            "domain": "continent_global", "model": model,
            "status": "not_estimable", "message": message,
        })

    write_tsv(
        pd.DataFrame(null_status),
        output_dir / f"{stem}_randomization_status.tsv",
        index=False,
    )

    summary = integrated_fungus_summary(
        selected, occurrences, host_rows, continent_rows, args.analysis_level
    )
    write_tsv(
        summary,
        output_dir / f"{rank_prefix}_selected_fungi_integrated_summary.tsv",
        index=False,
    )
    if args.analysis_level == "sh":
        extract_selected_fasta(
            selected_ids, confidence_fasta_files,
            output_dir / f"{rank_prefix}_selected_SH_sequences.fasta",
            output_dir / f"{rank_prefix}_selected_SH_fasta_status.tsv",
            str(input_confidence),
        )
        plot_selected_spatial_distributions(
            occurrences, selected_ids, args.world_geojson,
            output_dir / f"{rank_prefix}_selected_SH_spatial_maps",
        )
    append_log(analysis_log, "Completed input")


def _run_workflow_impl(analysis_level: str) -> int:
    args = parse_args(analysis_level)
    project_dir = Path(__file__).resolve().parent
    set_figure_font_scale(args.figure_font_scale)
    args.output_dir = args.output_dir.resolve()
    copy_globalfungi_inputs(args, project_dir)
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.r_script = args.r_script.resolve()
    args.world_geojson = (
        args.world_geojson
        or (Path(__file__).resolve().parent / "assets" / "ne_110m_admin_0_countries.geojson")
    ).resolve()
    if analysis_level == "sh" and not args.world_geojson.exists():
        raise FileNotFoundError(f"World boundary GeoJSON does not exist: {args.world_geojson}")
    require_r_environment(args)
    taxonomy_path = (args.taxonomy or (args.data_dir / "plant_taxonomy.tsv")).resolve()
    override_path = (
        args.taxonomy_overrides or (args.data_dir / "plant_taxonomy_overrides.tsv")
    ).resolve()
    paths = input_paths(args.data_dir, args.input_glob, [taxonomy_path, override_path], args.only)
    plant_names: list[str] = []
    valid_paths: list[Path] = []
    for path in paths:
        try:
            mapping = occurrence_column_map(path)
            host_column = mapping["host_candidate"]
            hosts = pd.read_csv(
                path, sep="\t", usecols=[host_column], encoding="utf-8-sig",
                low_memory=False,
            )[host_column]
        except Exception as exc:
            print(f"Skipping invalid unique-occurrence TSV {path.name}: {exc}", file=sys.stderr)
            continue
        valid_paths.append(path)
        plant_names.extend(
            genus for genus in hosts.map(host_plant_genus).unique() if not is_unknown(genus)
        )
    if not valid_paths:
        raise RuntimeError("No valid medium-or-higher unique-occurrence TSV files were found.")
    cpu_count = os.cpu_count() or 1
    args.effective_jobs = args.jobs if args.jobs > 0 else cpu_count
    args.effective_jobs = max(1, args.effective_jobs)
    args.effective_input_jobs = min(len(valid_paths), args.effective_jobs)
    write_run_configuration(args, valid_paths)
    write_tsv(
        pd.DataFrame(runtime_version_records(args)),
        args.output_dir / "runtime_versions.tsv",
        index=False,
    )
    taxonomy = prepare_taxonomy(
        plant_names, taxonomy_path, args.offline, args.taxonomy_workers,
        refresh=args.refresh_taxonomy, override_path=override_path,
    )
    if args.prepare_taxonomy_only:
        print(f"Wrote taxonomy cache: {taxonomy_path}")
        return 0
    fasta_files = fasta_paths(args.data_dir, args.fasta)
    if fasta_files:
        print("FASTA source(s): " + ", ".join(path.name for path in fasta_files))
    print(
        f"Parallel budget: {args.effective_jobs}; input workers: "
        f"{args.effective_input_jobs} "
        f"(inputs={len(valid_paths)}, CPUs={cpu_count}, no non-IQ-TREE cap)"
    )
    append_log(
        args.output_dir / "workflow.log",
        f"Parallel budget={args.effective_jobs}, input_workers={args.effective_input_jobs}",
    )
    if args.effective_input_jobs == 1:
        for path in valid_paths:
            process_one(path, taxonomy, fasta_files, args)
    else:
        try:
            with ProcessPoolExecutor(max_workers=args.effective_input_jobs) as pool:
                futures = {
                    pool.submit(process_one, path, taxonomy, fasta_files, args): path
                    for path in valid_paths
                }
                for future in as_completed(futures):
                    path = futures[future]
                    try:
                        future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            f"Parallel analysis failed for {path.name}: {exc}"
                        ) from exc
        except PermissionError as error:
            warnings.warn(
                "Process creation is unavailable; processing confidence inputs "
                f"sequentially ({error})."
            )
            append_log(
                args.output_dir / "workflow.log",
                f"Input-process fallback to sequential mode: {error}",
            )
            args.effective_input_jobs = 1
            write_tsv(
                pd.DataFrame([{
                    "reason": str(error),
                    "actual_input_workers": 1,
                }]),
                args.output_dir / "parallel_fallback.tsv",
                index=False,
            )
            for path in valid_paths:
                process_one(path, taxonomy, fasta_files, args)
    print(f"\nCompleted {len(valid_paths)} unique-occurrence file(s). Outputs: {args.output_dir}")
    append_log(args.output_dir / "workflow.log", f"Completed {len(valid_paths)} inputs")
    return 0


def run_workflow(analysis_level: str) -> int:
    """Time the complete workflow and preserve timing on success or failure."""
    RUN_STATE.clear()
    RUN_STATE.update({
        "started_at": datetime.now(timezone.utc).astimezone().isoformat(),
        "started_monotonic": time.perf_counter(),
        "finalized": False,
    })
    exit_code = 1
    try:
        exit_code = _run_workflow_impl(analysis_level)
    except BaseException:
        finalize_run_timing("failed", exit_code)
        raise
    else:
        finalize_run_timing(
            "completed" if exit_code == 0 else "nonzero_exit", exit_code
        )
        return exit_code


if __name__ == "__main__":
    raise SystemExit(
        "Run specificity_sh_workflow.py or specificity_genus_workflow.py, "
        "not specificity_common.py."
    )
