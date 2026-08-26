#!/usr/bin/env python3
"""Explore Helotiales-associated occurrence data downloaded from GlobalFungi.

Edit the USER SETTINGS below, or supply paths with command-line options.  The
workflow accepts CSV/TSV (including gzip-compressed) tables, maps common column
name variants to a standard schema, selects records for Helotiales or any other
named taxon, merges sample metadata, writes summaries, makes exploratory figures,
and optionally exports sequences.

Two occurrence layouts are supported:
1. A long table with one sample-SH/variant occurrence per row.
2. GlobalFungi's wide sample-by-SH abundance matrix. For this layout, supply
   either an SH taxonomy table (``--taxonomy``) or an explicit SH-ID list
   (``--sh-list``). Only selected SH columns are loaded, in chunks.

For every run, an additional exact ``sample_type == root`` subset is produced.
It includes reproducible host_candidate, host_confidence, associated_plants,
host_evidence, and source_study fields. If a FASTA or UNITE .tgz archive is
available, sample-SH FASTA records and matching tables are exported for all,
ITS1, ITS2, and ITSboth subsets.

Example for the GlobalFungi wide matrix::

    python helotiales_globalfungi_workflow.py \
      --occurrences GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz \
      --metadata GlobalFungi_5_sample_metadata.txt.gz \
      --taxonomy unite_sh_taxonomy.tsv \
      --target-taxon Helotiales \
      --output results_helotiales

Repeat ``--target-taxon`` to select several taxa. Any taxonomic rank represented
in the taxonomy table can be used (for example Agaricales, Pezizaceae, Mollisia).

Required packages:
    python -m pip install pandas numpy matplotlib

Optional packages:
    python -m pip install biopython

Biopython is used only for an input FASTA file. Country outlines are drawn from
the bundled Natural Earth GeoJSON file and require no additional Python package.

Scientific caveats
------------------
* GlobalFungi integrates studies using different primers, extraction methods,
  sequencing platforms, and depths. Read counts are therefore not directly
  comparable among studies.
* Presence/absence and sample-level occurrence are generally safer than reads
  for global distribution analyses.
* Web-interface genus searches can miss unidentified Helotiales variants.
* A publication-level analysis should ideally define Helotiales from current
  UNITE Species Hypothesis taxonomy, not only from a hand-curated genus list.
* Host specificity cannot be inferred where host metadata are unavailable.
* Apparent host specificity must be interpreted with geography, biome, sample
  type, pH, climate, study ID, primer region, and other sampling covariates.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import hashlib
import io
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
import unicodedata
import warnings
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Iterable, Optional


PROGRAM_STARTED_AT = datetime.now(timezone.utc).astimezone()
PROGRAM_STARTED_MONOTONIC = time.perf_counter()
RUN_LOG_STATE: dict[str, object] = {}

# A shared writable cache prevents repeated font-cache warnings on systems where
# the user's default Matplotlib cache is read-only.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "globalfungi_matplotlib"))

import matplotlib

matplotlib.use("Agg")  # Safe on servers and in non-interactive workflows.
matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.text import Text
import numpy as np
import pandas as pd


# =============================================================================
# 1. USER SETTINGS -- edit these values for an interactive/scripted run
# =============================================================================

OCCURRENCE_TABLE = Path("globalfungi_occurrences.tsv")
SAMPLE_METADATA_TABLE = Path("globalfungi_sample_metadata.tsv")
FASTA_FILE: Optional[Path] = None
TAXONOMY_TABLE: Optional[Path] = None
SH_LIST_FILE: Optional[Path] = None
OUTPUT_DIRECTORY = Path("helotiales_globalfungi_results")
WORLD_BOUNDARIES_FILE = (
    Path(__file__).resolve().parent / "assets" / "ne_110m_admin_0_countries.geojson"
)
FIGURE_FONT_SCALE = 1.0


def set_figure_font_scale(scale: float) -> None:
    """Set the multiplier applied to every text element before figure export."""
    global FIGURE_FONT_SCALE
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("--figure-font-scale must be a finite number greater than 0.")
    FIGURE_FONT_SCALE = float(scale)


def apply_figure_font_scale(fig: plt.Figure) -> None:
    """Apply the configured multiplier once to every text artist in a figure."""
    if not getattr(fig, "_globalfungi_font_scale_applied", False):
        for text in fig.findobj(match=Text):
            text.set_fontsize(text.get_fontsize() * FIGURE_FONT_SCALE)
        fig._globalfungi_font_scale_applied = True


def save_figure(fig: plt.Figure, output_path: Path, dpi: int = 300) -> None:
    """Save each figure as a vector PDF."""
    apply_figure_font_scale(fig)
    vector_path = output_path.with_suffix(".pdf")
    fig.savefig(vector_path, format="pdf", bbox_inches="tight")

# Taxon names can be orders, families, genera, or other names present in the
# taxonomy table. For example: ["Helotiales"], ["Agaricales"], or
# ["Helotiales", "Pezizales"]. Matching is case-insensitive.
TARGET_TAXA = ["Helotiales"]

HELOTIALES_GENERA = [
    "Oidiodendron",
    "Hyaloscypha",
    "Meliniomyces",
    "Phialocephala",
    "Cadophora",
    "Rhizoscyphus",
    "Lachnum",
    "Tetracladium",
    "Mollisia",
    "Leohumicola",
    "Cudoniella",
    "Vibrissea",
]

MINIMUM_READ_COUNT = 0
USE_PRESENCE_ABSENCE = True
WIDE_MATRIX_CHUNK_SIZE = 500
MAX_PARALLEL_JOBS = None
ROOT_SAMPLE_TYPES = ["root"]


# Aliases are compared after removing spaces/punctuation and ignoring case.
# Add download-specific names here without changing the rest of the workflow.
COLUMN_ALIASES = {
    "sample_id": [
        "sample_id", "SampleID", "sample", "sampleID", "sample_name",
        "sample accession", "sample_accession",
    ],
    "sequence_id": [
        "sequence_id", "variant_id", "sv_id", "molecular_taxon_id", "OTU",
        "OTU_ID", "ASV", "ASV_ID", "sequence variant", "feature_id",
    ],
    "reads": ["reads", "read_count", "read count", "abundance", "count", "n_reads"],
    "genus": ["genus", "Genus", "tax_genus"],
    "species": ["species", "Species", "tax_species"],
    "family": ["family", "Family", "tax_family"],
    "order": ["order", "Order", "tax_order"],
    "class": ["class", "Class", "tax_class"],
    "phylum": ["phylum", "Phylum", "tax_phylum"],
    "kingdom": ["kingdom", "Kingdom", "tax_kingdom"],
    "taxonomy": [
        "taxonomy", "Taxonomy", "taxonomic_annotation", "taxon",
        "lineage", "classification", "scientific_name",
    ],
    "sh_id": [
        "sh_id", "SH", "UNITE_SH", "UNITE_SH_ID", "species_hypothesis",
        "species hypothesis id", "unite species hypothesis",
    ],
    "latitude": ["latitude", "lat", "Latitude", "decimal_latitude", "decimalLatitude"],
    "longitude": [
        "longitude", "lon", "long", "Longitude", "decimal_longitude",
        "decimalLongitude", "lng",
    ],
    "country": ["country", "Country", "country_name"],
    "continent": ["continent", "Continent"],
    "biome": [
        "biome", "Biome", "biome_name", "ecosystem_classification",
        "environment_type",
    ],
    "sample_type": [
        "sample_type", "sample type", "SampleType", "material", "sample_material",
    ],
    "sample_type_detail": [
        "sample_type_detail", "sample_type_specification", "sample type specification",
    ],
    "substrate": ["substrate", "Substrate", "substrate_type"],
    "host_or_vegetation": [
        "host_or_vegetation", "host", "host_plant", "host plant", "plant_host",
        "vegetation", "vegetation_type", "dominant_vegetation",
    ],
    "dominant_plant_species": [
        "dominant_plant_species", "dominant plants", "dominant plant species",
    ],
    "associated_plants": ["associated_plants", "other_plant_species", "other plants"],
    "sampling_info": ["sampling_info", "sampling information"],
    "sample_description": [
        "sample_description", "sample_info", "sample description",
    ],
    "study_title": ["study_title", "paper_title", "title"],
    "study_doi": ["study_doi", "paper_doi", "doi"],
    "its_region": [
        "its_region", "ITS region", "marker", "target_region", "primer_region",
        "barcoding_region",
    ],
    "sequencing_platform": [
        "sequencing_platform", "sequencing platform", "platform", "instrument",
        "sequencing_instrument", "sequencing technology", "sequencing_technology",
    ],
    "sequence": ["sequence", "dna_sequence", "nucleotide_sequence", "DNA sequence"],
    "ph": ["ph", "pH", "soil_ph", "soil pH"],
    "mat": [
        "mat", "MAT", "mean_annual_temperature", "mean annual temperature",
        "annual_mean_temperature", "bio1", "MAT_study",
    ],
    "map": [
        "map", "MAP", "mean_annual_precipitation", "mean annual precipitation",
        "annual_precipitation", "bio12", "MAP_study",
    ],
    "soil_carbon": [
        "soil_carbon", "soil carbon", "soil_organic_carbon", "organic_carbon", "soc",
        "organic_C_content",
    ],
}

STANDARD_COLUMNS = list(COLUMN_ALIASES)
METADATA_PRIORITY_COLUMNS = {
    "latitude", "longitude", "country", "continent", "biome", "sample_type",
    "sample_type_detail", "substrate", "host_or_vegetation", "dominant_plant_species",
    "associated_plants", "sampling_info",
    "sample_description", "study_title", "study_doi", "ph", "mat", "map",
    "soil_carbon", "sequencing_platform",
}


def normalized_name(value: object) -> str:
    """Normalize a column name for conservative, case-insensitive matching."""
    return re.sub(r"[^a-z0-9]+", "", str(value).casefold())


def find_column(columns: Iterable[object], aliases: Iterable[str]) -> Optional[str]:
    """Return the first exact normalized alias match, or None."""
    lookup = {}
    for col in columns:
        lookup.setdefault(normalized_name(col), str(col))
    for alias in aliases:
        match = lookup.get(normalized_name(alias))
        if match is not None:
            return match
    return None


def detect_columns(df: pd.DataFrame) -> dict[str, Optional[str]]:
    """Map every standard field to the most likely source column."""
    return {
        standard: find_column(df.columns, aliases)
        for standard, aliases in COLUMN_ALIASES.items()
    }


def open_text(path: Path):
    """Open plain, gzip, or bzip2 text using UTF-8 with BOM handling."""
    if str(path).casefold().endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8-sig", newline="")
    if str(path).casefold().endswith((".bz2", ".bzip2")):
        return bz2.open(path, mode="rt", encoding="utf-8-sig", newline="")
    return path.open(mode="r", encoding="utf-8-sig", newline="")


def read_header(path: Path) -> tuple[list[str], str]:
    """Read only a table header and detect comma, tab, or semicolon delimiter."""
    if not path.exists():
        raise FileNotFoundError(f"Input table does not exist: {path}")
    with open_text(path) as handle:
        first_line = handle.readline()
    if not first_line:
        raise ValueError(f"Input table is empty: {path}")
    counts = {delimiter: first_line.count(delimiter) for delimiter in ["\t", ",", ";"]}
    delimiter = max(counts, key=counts.get)
    if counts[delimiter] == 0:
        raise ValueError(f"Could not detect a supported delimiter in: {path}")
    columns = next(csv.reader([first_line], delimiter=delimiter))
    return columns, delimiter


def load_table(path: Path, **kwargs) -> pd.DataFrame:
    """Load CSV/TSV text, including gzip files, with a fast explicit delimiter."""
    _, delimiter = read_header(path)
    return pd.read_csv(path, sep=delimiter, compression="infer", low_memory=False, **kwargs)


def load_taxonomy_table(path: Path) -> pd.DataFrame:
    """Load ordinary taxonomy tables or official headerless UNITE SHs.tax files."""
    columns, delimiter = read_header(path)
    headerless_unite = bool(columns) and bool(
        re.fullmatch(r"SH\d+(?:\.\d+FU)?", columns[0], flags=re.I)
    )
    if not headerless_unite:
        return load_table(path)

    # Official *.SHs.tax files have no header and contain:
    # SH, taxon node ID, kingdom, phylum, class, order, family, genus, species.
    names = [
        "sh_id", "taxon_id", "kingdom", "phylum", "class", "order",
        "family", "genus", "species",
    ]
    raw = pd.read_csv(
        path,
        sep=delimiter,
        compression="infer",
        header=None,
        names=names,
        usecols=range(min(len(columns), len(names))),
        low_memory=False,
    )
    rank_prefixes = [
        ("kingdom", "k__"), ("phylum", "p__"), ("class", "c__"),
        ("order", "o__"), ("family", "f__"), ("genus", "g__"),
        ("species", "s__"),
    ]
    raw["taxonomy"] = raw.apply(
        lambda row: ";".join(
            prefix + str(row[rank])
            for rank, prefix in rank_prefixes
            if rank in row.index and pd.notna(row[rank]) and str(row[rank]).strip()
        ),
        axis=1,
    )
    print("Detected official headerless UNITE SHs.tax taxonomy format.")
    return raw


def detect_occurrence_layout(path: Path) -> tuple[str, list[str], str]:
    """Recognize a long occurrence table or GlobalFungi's wide SH matrix."""
    columns, delimiter = read_header(path)
    if len(columns) > 2:
        candidate_columns = columns[1 : min(len(columns), 1001)]
        sh_fraction = np.mean(
            [bool(re.fullmatch(r"SH\d+(?:\.\d+FU)?", value, re.I)) for value in candidate_columns]
        )
        if sh_fraction >= 0.8:
            return "wide_sh_matrix", columns, delimiter
    return "long_occurrence", columns, delimiter


def inspect_table(name: str, df: pd.DataFrame, mapping: dict[str, Optional[str]]) -> None:
    """Print dimensions, columns, a preview, and detected standard fields."""
    print(f"\n{name}: {df.shape[0]:,} rows x {df.shape[1]:,} columns")
    print("Columns:", list(df.columns))
    print("Detected mapping:")
    for standard, source in mapping.items():
        if source is not None:
            print(f"  {standard:20s} <- {source}")
    print("First rows:")
    with pd.option_context("display.max_columns", 30, "display.width", 180):
        print(df.head())


def standardize_table(
    df: pd.DataFrame, mapping: dict[str, Optional[str]]
) -> pd.DataFrame:
    """Add standard columns while retaining every original input column."""
    out = df.copy()
    for standard in STANDARD_COLUMNS:
        source = mapping.get(standard)
        if source is None:
            out[standard] = pd.NA
        elif source != standard:
            out[standard] = out[source]

    # IDs are strings: preserving leading zeros is more useful than numeric IDs.
    for col in ["sample_id", "sequence_id", "sh_id"]:
        out[col] = out[col].astype("string").str.strip().replace("", pd.NA)
    for col in ["reads", "latitude", "longitude", "ph", "mat", "map", "soil_carbon"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    for col in [
        "genus", "species", "family", "order", "class", "phylum", "kingdom",
        "taxonomy", "country", "continent",
        "biome", "sample_type", "sample_type_detail", "substrate", "host_or_vegetation",
        "dominant_plant_species",
        "associated_plants", "sampling_info", "sample_description", "study_title",
        "study_doi", "its_region", "sequencing_platform",
        "sequence",
    ]:
        out[col] = out[col].astype("string").str.strip().replace("", pd.NA)
    return out


def validate_required_keys(
    occurrence_mapping: dict[str, Optional[str]], metadata_mapping: dict[str, Optional[str]]
) -> None:
    if occurrence_mapping["sample_id"] is None:
        raise ValueError(
            "No sample-ID column was detected in the occurrence table. Add its name "
            "to COLUMN_ALIASES['sample_id']."
        )
    if metadata_mapping["sample_id"] is None:
        raise ValueError(
            "No sample-ID column was detected in the metadata table. Add its name "
            "to COLUMN_ALIASES['sample_id']."
        )


def filter_minimum_reads(
    occurrence: pd.DataFrame,
    source_reads_column: Optional[str],
    minimum_reads: float,
) -> pd.DataFrame:
    """Apply a reads threshold only when a genuine reads column exists."""
    if source_reads_column is None:
        if minimum_reads > 0:
            warnings.warn("No reads column detected; minimum-read threshold was not applied.")
        return occurrence.copy()
    keep = occurrence["reads"].ge(minimum_reads).fillna(False)
    print(
        f"Read threshold >= {minimum_reads:g}: retained {int(keep.sum()):,} of "
        f"{len(occurrence):,} occurrence rows"
    )
    return occurrence.loc[keep].copy()


TAXONOMIC_RANK_COLUMNS = [
    "kingdom", "phylum", "class", "order", "family", "genus", "species"
]


def match_target_taxa(
    table: pd.DataFrame,
    targets: list[str],
    include_helotiales_genera: bool = True,
) -> tuple[pd.Series, pd.Series, pd.DataFrame]:
    """Match arbitrary taxon names against ranks and a taxonomy string.

    Returns a selection mask, semicolon-separated per-row reasons, and a
    criterion summary. Exact matches are used in rank columns; taxonomy-string
    matches are token bounded and case-insensitive.
    """
    criteria: dict[str, pd.Series] = {}
    index = table.index
    usable_targets = [str(target).strip() for target in targets if str(target).strip()]
    if not usable_targets:
        raise ValueError("At least one non-empty target taxon must be supplied.")

    for target in usable_targets:
        target_folded = target.casefold()
        for rank in TAXONOMIC_RANK_COLUMNS:
            if rank in table:
                criteria[f"{rank}_match:{target}"] = (
                    table[rank].astype("string").str.strip().str.casefold().eq(target_folded).fillna(False)
                )
        if "taxonomy" in table:
            criteria[f"taxonomy_match:{target}"] = table["taxonomy"].astype("string").str.contains(
                rf"(?<![A-Za-z0-9]){re.escape(target)}(?![A-Za-z0-9])",
                case=False,
                na=False,
                regex=True,
            )

    # Preserve the original exploratory genus-list fallback for Helotiales.
    if include_helotiales_genera and any(t.casefold() == "helotiales" for t in usable_targets):
        genus_set = {genus.casefold() for genus in HELOTIALES_GENERA}
        if "genus" in table:
            criteria["helotiales_genus_list_match"] = (
                table["genus"].astype("string").str.casefold().isin(genus_set)
            )
        if "taxonomy" in table:
            genus_pattern = "|".join(
                re.escape(genus) for genus in sorted(HELOTIALES_GENERA, key=len, reverse=True)
            )
            criteria["taxonomy_contains_helotiales_genus"] = (
                table["taxonomy"].astype("string").str.contains(
                    rf"(?<![A-Za-z])(?:{genus_pattern})(?![A-Za-z])",
                    case=False,
                    na=False,
                    regex=True,
                )
            )

    criterion_frame = pd.DataFrame(criteria, index=index).fillna(False)
    selected = criterion_frame.any(axis=1)
    reasons = pd.Series(pd.NA, index=index, dtype="string")
    reasons.loc[selected] = criterion_frame.loc[selected].apply(
        lambda row: ";".join(row.index[row.to_numpy(dtype=bool)]), axis=1
    )
    summary = pd.DataFrame(
        {
            "criterion": list(criteria),
            "n_records": [int(criterion_frame[name].sum()) for name in criteria],
        }
    )
    summary = pd.concat(
        [summary, pd.DataFrame([{
            "criterion": "selected_by_any_criterion",
            "n_records": int(selected.sum()),
        }])],
        ignore_index=True,
    )
    return selected, reasons, summary


def sh_base(value: object) -> Optional[str]:
    """Return version-independent UNITE SH stem (for example SH1234567)."""
    if pd.isna(value):
        return None
    match = re.match(r"^(SH\d+)", str(value).strip(), re.I)
    return match.group(1).upper() if match else None


def load_explicit_sh_list(path: Path) -> pd.DataFrame:
    """Read SH IDs from a one-per-line or delimited text file."""
    if not path.exists():
        raise FileNotFoundError(f"SH-list file does not exist: {path}")
    found: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            found.extend(re.findall(r"SH\d+(?:\.\d+FU)?", line, flags=re.I))
    unique = list(dict.fromkeys(value.upper() for value in found))
    if not unique:
        raise ValueError(f"No UNITE SH IDs were found in: {path}")
    return pd.DataFrame({
        "sh_id": unique,
        "target_taxon_match_reason": "explicit_sh_list",
    })


def select_shs_from_taxonomy(
    taxonomy_path: Path, targets: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load an SH taxonomy table and retain SHs assigned to target taxa."""
    raw = load_taxonomy_table(taxonomy_path)
    mapping = detect_columns(raw)
    inspect_table("SH taxonomy table", raw, mapping)
    if mapping["sh_id"] is None:
        raise ValueError(
            "No SH-ID column was detected in the taxonomy table. Add its column name "
            "to COLUMN_ALIASES['sh_id']."
        )
    taxonomy = standardize_table(raw, mapping)
    selected, reasons, summary = match_target_taxa(taxonomy, targets)
    result = taxonomy.loc[selected, STANDARD_COLUMNS].copy()
    result["target_taxon_match_reason"] = reasons.loc[selected]
    result = result[result["sh_id"].notna()].drop_duplicates("sh_id", keep="first")
    if result.empty:
        raise ValueError(
            f"No SH IDs matched target taxa {targets!r}. Check taxonomy columns and spelling."
        )
    return result, summary


def annotate_explicit_shs(
    selection: pd.DataFrame, taxonomy_path: Optional[Path]
) -> pd.DataFrame:
    """Optionally add taxonomy fields to explicitly selected SH IDs."""
    if taxonomy_path is None:
        for col in STANDARD_COLUMNS:
            if col not in selection:
                selection[col] = pd.NA
        return selection
    raw = load_taxonomy_table(taxonomy_path)
    mapping = detect_columns(raw)
    if mapping["sh_id"] is None:
        raise ValueError("No SH-ID column was detected in the taxonomy table.")
    taxonomy = standardize_table(raw, mapping)
    taxonomy["_sh_base"] = taxonomy["sh_id"].map(sh_base)
    selection["_sh_base"] = selection["sh_id"].map(sh_base)
    annotation_cols = ["_sh_base"] + [c for c in STANDARD_COLUMNS if c != "sh_id"]
    annotations = taxonomy[annotation_cols].drop_duplicates("_sh_base", keep="first")
    result = selection.merge(annotations, on="_sh_base", how="left")
    result.drop(columns="_sh_base", inplace=True)
    return result


def resolve_matrix_sh_columns(
    matrix_columns: list[str], selection: pd.DataFrame
) -> pd.DataFrame:
    """Resolve taxonomy SH IDs to matrix columns, tolerating version suffix changes."""
    sh_columns = [
        col for col in matrix_columns
        if re.fullmatch(r"SH\d+(?:\.\d+FU)?", col, flags=re.I)
    ]
    exact = {col.casefold(): col for col in sh_columns}
    by_base: dict[str, list[str]] = {}
    for col in sh_columns:
        base = sh_base(col)
        if base:
            by_base.setdefault(base, []).append(col)

    rows = []
    unresolved = 0
    for _, row in selection.iterrows():
        requested = str(row["sh_id"])
        matrix_id = exact.get(requested.casefold())
        if matrix_id is None:
            candidates = by_base.get(sh_base(requested) or "", [])
            if len(candidates) == 1:
                matrix_id = candidates[0]
        if matrix_id is None:
            unresolved += 1
            continue
        record = row.to_dict()
        record["taxonomy_sh_id"] = requested
        record["matrix_sh_id"] = matrix_id
        record["sh_id"] = matrix_id
        record["sequence_id"] = matrix_id
        rows.append(record)

    resolved = pd.DataFrame(rows).drop_duplicates("matrix_sh_id", keep="first")
    if unresolved:
        warnings.warn(
            f"{unresolved:,} selected taxonomy/list SH IDs were absent from the abundance matrix."
        )
    if resolved.empty:
        raise ValueError("None of the selected SH IDs occur as columns in the abundance matrix.")
    return resolved


def extract_wide_sh_occurrences(
    path: Path,
    matrix_columns: list[str],
    delimiter: str,
    resolved_taxonomy: pd.DataFrame,
    minimum_reads: float,
    chunk_size: int,
) -> pd.DataFrame:
    """Read only selected SH columns and convert non-zero cells to long records."""
    sample_source = find_column(matrix_columns, COLUMN_ALIASES["sample_id"])
    if sample_source is None:
        raise ValueError("No sample-ID column was detected in the wide abundance matrix.")
    selected_shs = resolved_taxonomy["matrix_sh_id"].astype(str).tolist()
    selected_set = set(selected_shs)
    usecols = lambda col: col == sample_source or col in selected_set
    chunks = pd.read_csv(
        path,
        sep=delimiter,
        compression="infer",
        usecols=usecols,
        chunksize=chunk_size,
        low_memory=False,
    )
    long_chunks: list[pd.DataFrame] = []
    processed = 0
    cutoff = max(float(minimum_reads), 0.0)
    print(
        f"\nWide SH matrix: reading {len(selected_shs):,} selected SH columns "
        f"from {len(matrix_columns) - 1:,} total SH columns"
    )
    for chunk in chunks:
        # Missing selected columns cannot occur after resolution, but reindex keeps
        # column order deterministic across pandas versions.
        values_frame = chunk.reindex(columns=selected_shs).apply(pd.to_numeric, errors="coerce")
        values = values_frame.to_numpy(dtype=float, na_value=0.0)
        keep = (values > 0) & (values >= cutoff)
        row_index, col_index = np.nonzero(keep)
        if len(row_index):
            long_chunks.append(pd.DataFrame({
                "sample_id": chunk[sample_source].astype("string").to_numpy()[row_index],
                "sequence_id": np.asarray(selected_shs, dtype=object)[col_index],
                "reads": values[row_index, col_index],
            }))
        processed += len(chunk)
        if processed % (chunk_size * 20) == 0:
            print(f"  processed {processed:,} samples", flush=True)

    if not long_chunks:
        raise ValueError(
            "Selected SH columns contained no abundance values above the requested threshold."
        )
    occurrence = pd.concat(long_chunks, ignore_index=True)
    annotation_columns = [
        "sequence_id", "sh_id", "taxonomy_sh_id", "taxonomy", "kingdom",
        "phylum", "class", "order", "family", "genus", "species",
        "its_region", "sequence", "target_taxon_match_reason",
    ]
    annotation = resolved_taxonomy[
        [col for col in annotation_columns if col in resolved_taxonomy]
    ].copy()
    occurrence = occurrence.merge(annotation, on="sequence_id", how="left", validate="many_to_one")
    occurrence["helotiales_match_reason"] = occurrence["target_taxon_match_reason"]
    print(
        f"Converted matrix to {len(occurrence):,} non-zero sample-SH occurrence records."
    )
    return occurrence


def extract_helotiales(
    occurrence: pd.DataFrame,
    genera: list[str],
    mapping: dict[str, Optional[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select rows by order, taxonomy text, or an editable genus list."""
    index = occurrence.index
    criteria: dict[str, pd.Series] = {}

    if mapping["order"] is not None:
        criteria["order_column_match"] = occurrence["order"].str.fullmatch(
            r"\s*Helotiales\s*", case=False, na=False
        )
    else:
        criteria["order_column_match"] = pd.Series(False, index=index)

    if mapping["taxonomy"] is not None:
        criteria["taxonomy_string_match"] = occurrence["taxonomy"].str.contains(
            r"\bHelotiales\b", case=False, na=False, regex=True
        )
    else:
        criteria["taxonomy_string_match"] = pd.Series(False, index=index)

    genus_set = {g.casefold() for g in genera}
    if mapping["genus"] is not None:
        criteria["genus_list_match"] = occurrence["genus"].str.casefold().isin(genus_set)
    else:
        criteria["genus_list_match"] = pd.Series(False, index=index)

    genus_pattern = "|".join(re.escape(g) for g in sorted(genera, key=len, reverse=True))
    if mapping["taxonomy"] is not None and genus_pattern:
        criteria["taxonomy_contains_helotiales_genus"] = occurrence["taxonomy"].str.contains(
            rf"(?<![A-Za-z])(?:{genus_pattern})(?![A-Za-z])",
            case=False,
            na=False,
            regex=True,
        )
    else:
        criteria["taxonomy_contains_helotiales_genus"] = pd.Series(False, index=index)

    criterion_frame = pd.DataFrame(criteria, index=index)
    selected = criterion_frame.any(axis=1)
    result = occurrence.loc[selected].copy()
    result["helotiales_match_reason"] = criterion_frame.loc[selected].apply(
        lambda row: ";".join(row.index[row.to_numpy(dtype=bool)]), axis=1
    )

    summary = pd.DataFrame(
        {
            "criterion": list(criteria),
            "n_records": [int(criteria[name].sum()) for name in criteria],
        }
    )
    summary = pd.concat(
        [
            summary,
            pd.DataFrame(
                [{"criterion": "selected_by_any_criterion", "n_records": int(selected.sum())}]
            ),
        ],
        ignore_index=True,
    )
    return result, summary


def deduplicate_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Reduce duplicate sample rows to one row, taking first non-null values."""
    valid = metadata[metadata["sample_id"].notna()].copy()
    n_duplicates = int(valid["sample_id"].duplicated(keep=False).sum())
    if n_duplicates:
        warnings.warn(
            f"Metadata contains {n_duplicates:,} rows belonging to duplicated sample IDs; "
            "first non-null value per column will be used for each sample."
        )
        valid = valid.groupby("sample_id", as_index=False, sort=False).first()
    return valid


def merge_occurrences_metadata(
    occurrence: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    """Left-join metadata, then coalesce standard columns by sensible priority."""
    metadata = deduplicate_metadata(metadata)
    merged = occurrence.merge(
        metadata,
        on="sample_id",
        how="left",
        suffixes=("__occ", "__meta"),
        validate="many_to_one",
        indicator="metadata_merge_status",
    )

    for col in STANDARD_COLUMNS:
        if col == "sample_id":
            continue
        occ_col, meta_col = f"{col}__occ", f"{col}__meta"
        if occ_col in merged and meta_col in merged:
            if col in METADATA_PRIORITY_COLUMNS:
                merged[col] = merged[meta_col].where(
                    merged[meta_col].notna(), merged[occ_col]
                )
            else:
                merged[col] = merged[occ_col].where(
                    merged[occ_col].notna(), merged[meta_col]
                )
            merged.drop(columns=[occ_col, meta_col], inplace=True)
        elif occ_col in merged:
            merged.rename(columns={occ_col: col}, inplace=True)
        elif meta_col in merged:
            merged.rename(columns={meta_col: col}, inplace=True)

    unmatched = int((merged["metadata_merge_status"] == "left_only").sum())
    if unmatched:
        warnings.warn(f"{unmatched:,} selected occurrence rows did not match sample metadata.")
    return merged


def valid_text(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().notna() & series.astype("string").str.strip().ne("")


def report_merged_data(df: pd.DataFrame, target_label: str = "Helotiales") -> pd.DataFrame:
    """Print and return core dataset diagnostics."""
    valid_coords = df["latitude"].between(-90, 90) & df["longitude"].between(-180, 180)
    host_or_substrate = valid_text(df["host_or_vegetation"]) | valid_text(df["substrate"])
    stats = {
        "occurrence_records": len(df),
        "unique_samples": df["sample_id"].nunique(dropna=True),
        "unique_sequence_variants": df["sequence_id"].nunique(dropna=True),
        "unique_unite_sh_ids": df["sh_id"].nunique(dropna=True),
        "countries": df["country"].nunique(dropna=True),
        "continents": df["continent"].nunique(dropna=True),
        "records_with_valid_coordinates": int(valid_coords.sum()),
        "records_with_host_or_substrate": int(host_or_substrate.sum()),
    }
    print(f"\nMerged target-taxon dataset ({target_label}):")
    for key, value in stats.items():
        print(f"  {key:36s} {value:,}")
    return pd.DataFrame({"metric": stats.keys(), "value": stats.values()})


def safe_nunique(series: pd.Series) -> int:
    return int(series.nunique(dropna=True))


def geographic_ranges_by_genus(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    work["genus"] = work["genus"].fillna("Unknown")
    work.loc[~work["latitude"].between(-90, 90), "latitude"] = np.nan
    work.loc[~work["longitude"].between(-180, 180), "longitude"] = np.nan
    return (
        work.groupby("genus", dropna=False)
        .agg(
            latitude_min=("latitude", "min"),
            latitude_max=("latitude", "max"),
            longitude_min=("longitude", "min"),
            longitude_max=("longitude", "max"),
            n_georeferenced_samples=(
                "sample_id",
                lambda s: s[work.loc[s.index, "latitude"].notna() & work.loc[s.index, "longitude"].notna()].nunique(),
            ),
        )
        .reset_index()
    )


def clean_metadata_text(value: object) -> Optional[str]:
    """Return useful metadata text, treating common missing-value tokens as absent."""
    if pd.isna(value):
        return None
    text = re.sub(r"\s+", " ", str(value)).strip()
    if not text or text.casefold() in {"na", "n/a", "none", "nan", "unknown", "not available"}:
        return None
    return text


def apply_host_metadata_priority(df: pd.DataFrame) -> pd.DataFrame:
    """Prefer direct host_or_vegetation; fall back to dominant_plant_species."""
    out = df.copy()
    direct = out["host_or_vegetation"].map(clean_metadata_text).astype("string")
    dominant = out["dominant_plant_species"].map(clean_metadata_text).astype("string")
    out["host_or_vegetation_direct"] = direct
    out["host_or_vegetation"] = direct.fillna(dominant)
    out["host_or_vegetation_source"] = pd.Series(
        np.select(
            [direct.notna(), dominant.notna()],
            ["host_or_vegetation", "dominant_plant_species_fallback"],
            default="unknown",
        ),
        index=out.index,
        dtype="string",
    )
    return out


def join_unique_text(values: pd.Series, separator: str = "; ") -> object:
    cleaned = [clean_metadata_text(value) for value in values]
    unique = list(dict.fromkeys(value for value in cleaned if value is not None))
    return separator.join(unique) if unique else pd.NA


def host_evidence_text(row: pd.Series) -> str:
    """Build a transparent evidence synopsis from GlobalFungi metadata fields."""
    parts = []
    for column, label in [
        ("host_or_vegetation_source", "host_candidate_source"),
        ("sample_type", "sample_type"),
        ("sample_type_detail", "sample_type_detail"),
        ("sample_description", "sample_description"),
        ("sampling_info", "sampling_info"),
    ]:
        value = clean_metadata_text(row.get(column))
        if value:
            parts.append(f"{label}: {value}")
    return " | ".join(parts)[:1000] if parts else "metadata evidence unavailable"


def source_study_text(row: pd.Series) -> object:
    title = clean_metadata_text(row.get("study_title"))
    doi = clean_metadata_text(row.get("study_doi"))
    if title and doi:
        return f"{title} | DOI: {doi}"
    return title or doi or pd.NA


def host_confidence_for_row(row: pd.Series) -> str:
    """Apply the requested reproducible high/medium/low/unknown rules."""
    host = clean_metadata_text(row.get("host_candidate"))
    if host is None:
        return "unknown"

    evidence = " ".join(
        filter(
            None,
            [
                clean_metadata_text(row.get("sample_type")),
                clean_metadata_text(row.get("sample_type_detail")),
                clean_metadata_text(row.get("sample_description")),
                clean_metadata_text(row.get("sampling_info")),
            ],
        )
    ).casefold()
    multiple_host = bool(
        re.search(r"[,;/|]", host)
        or re.search(r"\b(?:and|mixed|multiple|various|several)\b", host, flags=re.I)
    )
    root_explicit = bool(
        re.search(r"\b(?:root|roots|fine root|root tip|root tips|root fragment|root system)\b", evidence)
    )
    focal_or_separated = bool(
        re.search(
            r"(?:focal|target plant|host plant|individual plant|individual tree|each plant|"
            r"each tree|per tree|per each tree|root separation|separated roots|washed roots|"
            r"surface[- ]sterili[sz]ed|root tip samples|roots from)",
            evidence,
        )
    )
    vegetation_level_only = bool(
        re.search(r"\b(?:stand|vegetation|plant community|forest plot|mixed forest)\b", evidence)
        and not focal_or_separated
    )

    if multiple_host or vegetation_level_only:
        return "low"
    if root_explicit and focal_or_separated:
        return "high"
    return "medium"


def strongest_host_confidence(values: pd.Series) -> str:
    """Use the strongest evidence when duplicate location-host records disagree."""
    priority = {"unknown": 0, "low": 1, "medium": 2, "high": 3}
    cleaned = [str(value).strip().casefold() for value in values.dropna()]
    cleaned = [value if value in priority else "unknown" for value in cleaned]
    return max(cleaned or ["unknown"], key=priority.get)


def join_delimited_unique(values: pd.Series) -> object:
    """Flatten semicolon-delimited member lists and retain first-seen order."""
    items: list[str] = []
    for value in values.dropna():
        items.extend(part.strip() for part in str(value).split(";") if part.strip())
    unique = list(dict.fromkeys(items))
    return "; ".join(unique) if unique else pd.NA


def collapse_unique_occurrences(observations: pd.DataFrame) -> pd.DataFrame:
    """Collapse equal SH-latitude-longitude-host records into unique occurrences.

    Rows missing any of latitude, longitude, or host_candidate remain separate
    sample-SH observations because their equality cannot be established.
    """
    if observations.empty:
        return observations.copy()
    work = observations.copy()
    work["latitude"] = pd.to_numeric(work["latitude"], errors="coerce")
    work["longitude"] = pd.to_numeric(work["longitude"], errors="coerce")
    work["host_candidate"] = work["host_candidate"].map(clean_metadata_text).astype("string")
    work["_source_order"] = np.arange(len(work))
    complete = (
        work["latitude"].notna()
        & work["longitude"].notna()
        & work["host_candidate"].notna()
    )
    keys = ["sh_id", "latitude", "longitude", "host_candidate"]
    complete_work = work.loc[complete].copy()
    grouped = complete_work.groupby(keys, sort=False, dropna=False)
    collapsed = complete_work.drop_duplicates(keys, keep="first").set_index(keys)
    collapsed["host_confidence"] = grouped["host_confidence"].agg(
        strongest_host_confidence
    )
    collapsed["host_confidences"] = grouped["host_confidence"].agg(join_unique_text)
    collapsed["n_sample_sh_occurrences"] = grouped.size()
    collapsed["n_sample_ids"] = grouped["sample_id"].nunique(dropna=True)
    collapsed["sample_ids"] = grouped["sample_id"].agg(join_unique_text)
    collapsed["n_source_studies"] = grouped["source_study"].nunique(dropna=True)
    collapsed["source_studies"] = grouped["source_study"].agg(join_unique_text)
    collapsed["n_source_sequence_records"] = grouped[
        "n_source_sequence_records"
    ].sum(min_count=1)
    collapsed["member_sequence_ids"] = grouped["member_sequence_ids"].agg(
        join_delimited_unique
    )
    collapsed["n_member_sequence_ids"] = collapsed["member_sequence_ids"].map(
        lambda value: len(str(value).split(";")) if pd.notna(value) else 0
    )
    collapsed["member_its_regions"] = grouped["member_its_regions"].agg(
        join_delimited_unique
    )
    collapsed["member_sequencing_platforms"] = grouped[
        "member_sequencing_platforms"
    ].agg(join_delimited_unique)
    collapsed["total_reads_in_occurrence"] = grouped[
        "total_reads_in_occurrence"
    ].sum(min_count=1)
    collapsed["_source_order"] = grouped["_source_order"].min()
    collapsed["occurrence_count"] = 1
    collapsed["unique_occurrence_count"] = 1
    collapsed["unique_occurrence_basis"] = (
        "same_sh_latitude_longitude_host_candidate"
    )
    collapsed = collapsed.reset_index()

    incomplete = work.loc[~complete].copy()
    incomplete["host_confidences"] = incomplete["host_confidence"]
    incomplete["n_sample_sh_occurrences"] = 1
    incomplete["n_sample_ids"] = incomplete["sample_id"].notna().astype(int)
    incomplete["sample_ids"] = incomplete["sample_id"]
    incomplete["n_source_studies"] = incomplete["source_study"].notna().astype(int)
    incomplete["source_studies"] = incomplete["source_study"]
    incomplete["occurrence_count"] = 1
    incomplete["unique_occurrence_count"] = 1
    incomplete["unique_occurrence_basis"] = "incomplete_key_retained_by_sample"

    unique = pd.concat([collapsed, incomplete], ignore_index=True, sort=False)
    unique = unique.sort_values("_source_order", kind="stable").drop(
        columns="_source_order"
    )
    leading = [
        "sh_id", "latitude", "longitude", "host_candidate", "host_confidence",
        "unique_occurrence_count", "n_sample_sh_occurrences", "sample_id",
        "n_sample_ids", "sample_ids", "source_study", "n_source_studies",
        "source_studies", "host_confidences", "unique_occurrence_basis",
    ]
    return unique[
        leading + [column for column in unique.columns if column not in leading]
    ].reset_index(drop=True)


def build_root_host_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create annotated root records and unique SH-location-host occurrences.

    Sequence records first collapse to sample-SH observations. Observations with
    equal SH, latitude, longitude, and host_candidate then collapse to one unique
    occurrence, while incomplete location/host keys remain sample-specific.
    """
    root_types = {value.casefold() for value in ROOT_SAMPLE_TYPES}
    root_mask = df["sample_type"].astype("string").str.strip().str.casefold().isin(root_types)
    root = df.loc[root_mask].copy()
    if root.empty:
        return root, pd.DataFrame()

    root["host_candidate"] = root["host_or_vegetation"].map(clean_metadata_text).astype("string")
    root["host_candidate_source"] = root["host_or_vegetation_source"].astype("string")
    root["associated_plants"] = root["associated_plants"].map(clean_metadata_text).astype("string")
    root["host_evidence"] = root.apply(host_evidence_text, axis=1)
    root["source_study"] = root.apply(source_study_text, axis=1)
    root["host_confidence"] = root.apply(host_confidence_for_row, axis=1)

    work = root.copy()
    observation_keys = [
        "source_study", "sample_id", "sh_id", "host_candidate", "host_confidence"
    ]
    for column in observation_keys:
        work[column] = work[column].astype("string").fillna("Unknown")
    grouped = work.groupby(observation_keys, sort=False, dropna=False)
    observations = work.drop_duplicates(observation_keys, keep="first").copy()
    observations = observations.set_index(observation_keys)
    observations["occurrence_count"] = 1
    observations["n_source_sequence_records"] = grouped.size()
    observations["n_member_sequence_ids"] = grouped["sequence_id"].nunique(dropna=True)
    observations["member_sequence_ids"] = grouped["sequence_id"].agg(join_unique_text)
    observations["member_its_regions"] = grouped["its_region"].agg(join_unique_text)
    observations["member_sequencing_platforms"] = grouped[
        "sequencing_platform"
    ].agg(join_unique_text)
    observations["total_reads_in_occurrence"] = grouped["reads"].sum(min_count=1)
    observations = observations.reset_index()
    leading = [
        "sh_id", "sample_id", "source_study", "host_candidate", "host_confidence",
        "occurrence_count", "n_source_sequence_records", "n_member_sequence_ids",
        "member_sequence_ids", "member_its_regions", "member_sequencing_platforms",
        "total_reads_in_occurrence",
    ]
    observations = observations[
        leading + [column for column in observations.columns if column not in leading]
    ]
    return root, collapse_unique_occurrences(observations)


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    # Keep every logical record on exactly one physical line. Long free-text
    # metadata can contain tabs or Unicode/newline separators in some studies.
    clean = df.copy()
    text_columns = clean.select_dtypes(include=["object", "string"]).columns
    punctuation_translation = str.maketrans({
        "‐": "-", "‑": "-", "‒": "-", "–": "-", "—": "-", "―": "-",
        "−": "-", "“": '"', "”": '"', "„": '"', "‘": "'", "’": "'",
        "×": "x", " ": " ",
    })
    for column in text_columns:
        clean[column] = clean[column].map(
            lambda value: re.sub(
                r"[\t\r\n\u2028\u2029]+",
                " ",
                unicodedata.normalize("NFC", value).translate(punctuation_translation),
            )
            if isinstance(value, str)
            else value
        )
    # Write canonical tab-separated UTF-8 text without a BOM. Some software
    # treats a BOM-prefixed .tsv as generic text instead of a TSV file.
    # Converting punctuation to ASCII above prevents the known mojibake case.
    # mis-detection (for example, "Land‐use" appearing as "Land窶訊se").
    clean.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
        lineterminator="\n",
        encoding="utf-8",
    )
    print(f"Wrote {path}")


def initialize_run_log(
    output_dir: Path, args: argparse.Namespace, effective_jobs: int
) -> None:
    """Write the exact command and resolved run configuration."""
    command = shlex.join([sys.executable, *sys.argv])
    command_path = output_dir / "run_command.sh"
    with command_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        handle.write(command + "\n")
    settings = vars(args).copy()
    settings.update({
        "command": command,
        "parallel_jobs_resolved": effective_jobs,
        "occurrence_unit": "unique SH-latitude-longitude-host_candidate",
    })
    save_tsv(
        pd.DataFrame([
            {"setting": key, "value": str(value)}
            for key, value in sorted(settings.items())
        ]),
        output_dir / "run_configuration.tsv",
    )
    log_path = output_dir / "workflow.log"
    with log_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"start_time\t{PROGRAM_STARTED_AT.isoformat()}\n")
        handle.write(f"command\t{command}\n")
        handle.write(f"parallel_jobs_resolved\t{effective_jobs}\n")
        handle.write("occurrence_unit\tunique SH-latitude-longitude-host_candidate\n")
    RUN_LOG_STATE.update({
        "output_dir": output_dir,
        "log_path": log_path,
        "command": command,
        "finalized": False,
    })


def finalize_run_log(status: str, exit_code: int) -> None:
    """Record completion status and elapsed wall-clock calculation time."""
    if not RUN_LOG_STATE or RUN_LOG_STATE.get("finalized"):
        return
    ended_at = datetime.now(timezone.utc).astimezone()
    elapsed = time.perf_counter() - PROGRAM_STARTED_MONOTONIC
    log_path = Path(RUN_LOG_STATE["log_path"])
    with log_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"end_time\t{ended_at.isoformat()}\n")
        handle.write(f"elapsed_seconds\t{elapsed:.6f}\n")
        handle.write(f"status\t{status}\n")
        handle.write(f"exit_code\t{exit_code}\n")
    timing = pd.DataFrame([{
        "start_time": PROGRAM_STARTED_AT.isoformat(),
        "end_time": ended_at.isoformat(),
        "elapsed_seconds": elapsed,
        "status": status,
        "exit_code": exit_code,
        "command": RUN_LOG_STATE["command"],
    }])
    save_tsv(timing, Path(RUN_LOG_STATE["output_dir"]) / "run_timing.tsv")
    RUN_LOG_STATE["finalized"] = True


def resolve_parallel_jobs(requested: int) -> int:
    """Resolve --jobs; zero means all available logical CPUs."""
    if requested < 0:
        raise ValueError("--jobs must be 0 (automatic) or a positive integer.")
    if requested == 0:
        return max(1, os.cpu_count() or 1)
    return requested


def runtime_version_records(
    script_path: Path,
    requested_jobs: int,
    effective_jobs: int,
) -> list[dict[str, str]]:
    """Inventory the runtime, all installed Python/R packages, and run settings."""
    records: list[dict[str, str]] = []

    def add(
        category: str,
        name: str,
        version: object = "",
        path: object = "",
        details: object = "",
    ) -> None:
        records.append({
            "category": str(category),
            "name": str(name),
            "version": str(version),
            "path": str(path),
            "details": str(details),
        })

    add("run", "timestamp", datetime.now(timezone.utc).astimezone().isoformat())
    add("run", "command", details=" ".join([sys.executable, *sys.argv]))
    add("run", "parallel_jobs_requested", requested_jobs)
    add("run", "parallel_jobs_resolved", effective_jobs)
    add("system", "platform", platform.platform())
    add("system", "machine", platform.machine())
    add("system", "processor", platform.processor())
    add(
        "python",
        platform.python_implementation(),
        platform.python_version(),
        sys.executable,
        sys.version.replace("\n", " "),
    )
    if script_path.exists():
        digest = hashlib.sha256(script_path.read_bytes()).hexdigest()
        add("script", script_path.name, digest, script_path, "SHA-256")

    python_packages = []
    for distribution in importlib_metadata.distributions():
        name = distribution.metadata.get("Name") or "Unknown"
        python_packages.append((name, distribution.version, distribution.locate_file("")))
    for name, version, location in sorted(
        python_packages,
        key=lambda item: (str(item[0]).casefold(), str(item[1]), str(item[2])),
    ):
        add("python_package", name, version, location)

    rscript = shutil.which("Rscript")
    if rscript is None:
        add("R", "Rscript", "not_found", details="R is not installed or not on PATH")
        return records
    try:
        version_result = subprocess.run(
            [rscript, "--version"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        r_version = (version_result.stdout or version_result.stderr).strip()
        add("R", "R", r_version, rscript)
        package_result = subprocess.run(
            [
                rscript,
                "-e",
                "ip <- installed.packages()[,c('Package','Version','LibPath'),drop=FALSE]; "
                "write.table(ip, sep='\\t', row.names=FALSE, col.names=FALSE, quote=FALSE)",
            ],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        if package_result.returncode == 0:
            for line in package_result.stdout.splitlines():
                fields = line.split("\t", 2)
                if len(fields) >= 2:
                    add(
                        "R_package",
                        fields[0],
                        fields[1],
                        fields[2] if len(fields) == 3 else "",
                    )
        else:
            add("R", "installed_packages", "unavailable", details=package_result.stderr.strip())
    except (OSError, subprocess.SubprocessError) as error:
        add("R", "Rscript", "inspection_failed", rscript, error)
    return records


def validate_root_occurrences(df: pd.DataFrame) -> pd.DataFrame:
    """Validate unique SH-location-host occurrences and representative sample IDs."""
    out = df.copy()
    sh_pattern = r"SH\d+(?:\.\d+FU)?"
    sh_valid = out["sh_id"].astype("string").str.fullmatch(sh_pattern, na=False)
    sequence_valid = out["sequence_id"].astype("string").str.fullmatch(
        sh_pattern, na=False
    )
    out.loc[~sh_valid & sequence_valid, "sh_id"] = out.loc[
        ~sh_valid & sequence_valid, "sequence_id"
    ]
    sh_valid = out["sh_id"].astype("string").str.fullmatch(sh_pattern, na=False)
    sample_valid = out["sample_id"].astype("string").str.strip().ne("").fillna(False)
    invalid = ~(sh_valid & sample_valid)
    if invalid.any():
        examples = out.loc[invalid, ["sh_id", "sample_id"]].head(5).to_dict("records")
        raise ValueError(
            f"Root occurrence table contains {int(invalid.sum()):,} row(s) without a "
            f"valid SH ID or sample_id; examples: {examples}"
        )
    complete = (
        pd.to_numeric(out["latitude"], errors="coerce").notna()
        & pd.to_numeric(out["longitude"], errors="coerce").notna()
        & out["host_candidate"].map(clean_metadata_text).notna()
    )
    duplicate_keys = ["sh_id", "latitude", "longitude", "host_candidate"]
    duplicated = complete & out.duplicated(duplicate_keys, keep=False)
    if duplicated.any():
        raise ValueError(
            f"Root occurrence table contains {int(duplicated.sum()):,} duplicated "
            "complete SH-latitude-longitude-host unique-occurrence keys."
        )
    first_columns = [
        "sh_id", "latitude", "longitude", "host_candidate", "host_confidence",
        "sample_id",
    ]
    return out[first_columns + [column for column in out.columns if column not in first_columns]]


def verify_tsv_data_lines_start_with_sh(path: Path) -> None:
    """Verify SH IDs, TSV field counts, and Excel-compatible cell lengths."""
    pattern = re.compile(r"^SH\d+(?:\.\d+FU)?\t")
    with path.open(mode="r", encoding="utf-8") as handle:
        header = next(handle, None)  # Header is intentionally not an SH record.
        if header is None:
            raise RuntimeError(f"TSV output is empty: {path}")
        expected_tabs = header.rstrip("\n").count("\t")
        if expected_tabs == 0:
            raise RuntimeError(f"TSV header contains no tab separators: {path}")
        for line_number, line in enumerate(handle, start=2):
            if not pattern.match(line):
                raise RuntimeError(
                    f"Invalid physical line {line_number} in {path}: it does not start with an SH ID."
                )
            actual_tabs = line.rstrip("\n").count("\t")
            if actual_tabs != expected_tabs:
                raise RuntimeError(
                    f"Invalid TSV field count on physical line {line_number} in {path}: "
                    f"expected {expected_tabs} tab separator(s), found {actual_tabs}."
                )
    with path.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header_fields = next(reader)
        for record_number, fields in enumerate(reader, start=2):
            if len(fields) != len(header_fields):
                raise RuntimeError(
                    f"Invalid parsed field count on record {record_number} in {path}: "
                    f"expected {len(header_fields)}, found {len(fields)}."
                )
            oversized = [header_fields[i] for i, value in enumerate(fields) if len(value) > 32_767]
            if oversized:
                raise RuntimeError(
                    f"Excel cell-length limit exceeded on record {record_number} in {path}: "
                    f"{', '.join(oversized)}."
                )


def draw_world_boundaries(ax: plt.Axes, geojson_path: Path) -> int:
    """Draw bundled Natural Earth country/coast lines without GeoPandas."""
    if not geojson_path.is_file():
        raise FileNotFoundError(
            f"World-boundary GeoJSON does not exist: {geojson_path}. "
            "Restore assets/ne_110m_admin_0_countries.geojson or pass --world-geojson."
        )
    with geojson_path.open(encoding="utf-8") as handle:
        features = json.load(handle).get("features", [])

    segments: list[np.ndarray] = []
    for feature in features:
        geometry = feature.get("geometry") or {}
        geometry_type = geometry.get("type")
        coordinates = geometry.get("coordinates") or []
        polygons = [coordinates] if geometry_type == "Polygon" else coordinates
        if geometry_type not in {"Polygon", "MultiPolygon"}:
            continue
        for polygon in polygons:
            for ring in polygon:
                points = np.asarray(ring, dtype=float)
                if points.ndim != 2 or points.shape[0] < 2 or points.shape[1] < 2:
                    continue
                # Do not draw a line across the whole map where a polygon crosses
                # the antimeridian.
                breaks = np.flatnonzero(np.abs(np.diff(points[:, 0])) > 180) + 1
                for part in np.split(points[:, :2], breaks):
                    if len(part) >= 2:
                        segments.append(part)
    if not segments:
        raise ValueError(f"No polygon boundaries were found in {geojson_path}.")
    ax.add_collection(
        LineCollection(
            segments,
            colors="#686868",
            linewidths=0.45,
            alpha=0.9,
            zorder=1,
        )
    )
    return len(segments)


def plot_world_points(
    df: pd.DataFrame,
    path: Path,
    world_geojson: Path,
    target_label: str = "Helotiales-associated",
    occurrence_count: Optional[int] = None,
    occurrence_label: str = "Occurrences",
) -> None:
    valid = df.loc[
        df["latitude"].between(-90, 90) & df["longitude"].between(-180, 180),
        ["latitude", "longitude"],
    ].copy()
    coords = valid.drop_duplicates()
    if coords.empty:
        print("Skipping world map: no valid coordinates.")
        return
    n_sites = len(coords)
    n_occurrences = int(len(df) if occurrence_count is None else occurrence_count)
    fig, ax = plt.subplots(figsize=(12, 6.2))
    ax.set_facecolor("#f7fbfc")
    n_boundary_segments = draw_world_boundaries(ax, world_geojson)
    ax.scatter(
        coords["longitude"], coords["latitude"], s=10, alpha=0.55,
        color="#8b2f67", edgecolors="none", zorder=2,
    )
    ax.set(xlim=(-180, 180), ylim=(-90, 90), xlabel="Longitude", ylabel="Latitude")
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.grid(color="#aaaaaa", linewidth=0.35, alpha=0.5)
    ax.set_title(target_label)
    ax.text(
        0.015,
        0.04,
        f"Sampling sites = {n_sites:,}\n{occurrence_label} = {n_occurrences:,}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#808080",
            "linewidth": 0.6,
            "alpha": 0.88,
        },
        zorder=3,
    )
    apply_figure_font_scale(fig)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)
    print(
        f"World map: drew {n_boundary_segments:,} land-boundary segments from "
        f"{world_geojson}."
    )


def plot_latitude_histogram(df: pd.DataFrame, path: Path) -> None:
    latitude = df.loc[df["latitude"].between(-90, 90), ["sample_id", "latitude"]].drop_duplicates()["latitude"]
    if latitude.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(latitude, bins=np.arange(-90, 95, 5), color="#516a9b", edgecolor="white")
    ax.set(xlabel="Latitude", ylabel="Number of unique samples", title="Latitude distribution")
    apply_figure_font_scale(fig)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def parse_unite_fasta_annotation(header: str) -> dict[str, str]:
    annotation: dict[str, str] = {}
    for rank, value in re.findall(r"(?:^|;)([kpcofgs])__([^;|]+)", header):
        column = {
            "k": "kingdom", "p": "phylum", "c": "class", "o": "order",
            "f": "family", "g": "genus", "s": "species",
        }[rank]
        annotation[column] = value.replace("_", " ")
    return annotation


SEQUENCING_PLATFORM_ORDER = {
    "Sanger": 0,
    "PacBio": 1,
    "Illumina": 2,
    "DNBSEQ": 3,
    "Ion Torrent": 4,
    "454 Roche": 5,
    "Other": 6,
    "Unknown": 7,
}


def normalized_sequencing_platform(value: object) -> str:
    """Map platform descriptions to the requested accuracy-priority groups."""
    text = clean_metadata_text(value)
    if text is None:
        return "Unknown"
    normalized = normalized_name(text)
    if any(token in normalized for token in ["sanger", "capillary", "abi3130", "abi3500", "abi3730"]):
        return "Sanger"
    if any(token in normalized for token in ["pacbio", "smrt", "sequel", "revio", "hifi"]):
        return "PacBio"
    if any(
        token in normalized
        for token in ["illumina", "miseq", "hiseq", "nextseq", "novaseq", "miniseq", "iseq"]
    ):
        return "Illumina"
    if any(token in normalized for token in ["dnbseq", "bgiseq", "mgiseq", "mgi"]):
        return "DNBSEQ"
    if "iontorrent" in normalized or "ionproton" in normalized or "ionpgm" in normalized:
        return "Ion Torrent"
    if "454" in normalized or "roche454" in normalized:
        return "454 Roche"
    return "Other"


def sequencing_platform_rank(value: object) -> int:
    return SEQUENCING_PLATFORM_ORDER[normalized_sequencing_platform(value)]


def sequence_length(sequence: str) -> int:
    """Count called bases, excluding alignment gaps and placeholder punctuation."""
    return len(re.sub(r"[-?.]", "", sequence))


def best_sequencing_platform(values: Iterable[object]) -> str:
    platforms = [normalized_sequencing_platform(value) for value in values]
    if not platforms:
        return "Unknown"
    return min(platforms, key=lambda value: (SEQUENCING_PLATFORM_ORDER[value], value))


def sequencing_platform_from_fasta_header(header: str) -> str:
    """Recognize an explicitly named platform without treating arbitrary text as Other."""
    platform_name = normalized_sequencing_platform(header)
    if platform_name != "Other":
        return platform_name
    normalized = normalized_name(header)
    if any(
        token in normalized
        for token in ["nanopore", "minion", "promethion"]
    ):
        return "Other"
    return "Unknown"


def sequence_selection_key(sequence: str, platform_name: object, identifier: object) -> tuple:
    """Longest sequence first; platform priority and identifier break ties."""
    return (
        -sequence_length(sequence),
        sequencing_platform_rank(platform_name),
        str(identifier).casefold(),
    )


def iter_fasta(handle):
    header: Optional[str] = None
    sequence_parts: list[str] = []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence_parts).upper()
            header = line[1:]
            sequence_parts = []
        else:
            sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts).upper()


def read_fasta(
    path: Path, target_ids: Optional[set[str]] = None
) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Read only target records from FASTA or a UNITE .tgz/.tar.gz archive."""
    if not path.exists():
        raise FileNotFoundError(f"FASTA file does not exist: {path}")
    target_ids = {str(value) for value in target_ids} if target_ids else set()
    target_by_base = {sh_base(value): value for value in target_ids if sh_base(value)}
    sequences: dict[str, str] = {}
    annotations: dict[str, dict[str, str]] = {}
    duplicate_candidates = 0

    def consume(handle) -> None:
        nonlocal duplicate_candidates
        for header, sequence in iter_fasta(handle):
            sh_match = re.search(r"SH\d+(?:\.\d+FU)?", header, flags=re.I)
            record_id = sh_match.group(0).upper() if sh_match else header.split()[0]
            output_id = record_id
            if target_ids:
                if record_id in target_ids:
                    output_id = record_id
                else:
                    output_id = target_by_base.get(sh_base(record_id))
                if output_id is None:
                    continue
            annotation = parse_unite_fasta_annotation(header)
            header_parts = header.split("|")
            annotation["accession"] = header_parts[1] if len(header_parts) > 1 else header.split()[0]
            annotation["sequencing_platform"] = sequencing_platform_from_fasta_header(header)
            annotation["fasta_header"] = header
            annotation["fasta_candidate_count"] = "1"
            if output_id in sequences:
                duplicate_candidates += 1
                previous = annotations[output_id]
                annotation["fasta_candidate_count"] = str(
                    int(previous.get("fasta_candidate_count", "1")) + 1
                )
                new_key = sequence_selection_key(
                    sequence, annotation["sequencing_platform"], annotation["accession"]
                )
                old_key = sequence_selection_key(
                    sequences[output_id],
                    previous.get("sequencing_platform", "Unknown"),
                    previous.get("accession", output_id),
                )
                if new_key >= old_key:
                    previous["fasta_candidate_count"] = annotation["fasta_candidate_count"]
                    continue
            sequences[output_id] = sequence
            annotations[output_id] = annotation

    lower_name = path.name.casefold()
    if lower_name.endswith((".tgz", ".tar.gz", ".tar")):
        with tarfile.open(path, mode="r:*") as archive:
            members = [
                member for member in archive.getmembers()
                if member.isfile() and member.name.casefold().endswith((".fa", ".fas", ".fasta"))
            ]
            preferred = [member for member in members if "_dev." not in member.name.casefold()]
            chosen = preferred[0] if preferred else (members[0] if members else None)
            if chosen is None:
                raise ValueError(f"No FASTA file was found inside archive: {path}")
            print(f"Reading UNITE FASTA archive member: {chosen.name}")
            binary = archive.extractfile(chosen)
            if binary is None:
                raise ValueError(f"Could not read FASTA member: {chosen.name}")
            with io.TextIOWrapper(binary, encoding="utf-8") as handle:
                consume(handle)
    elif lower_name.endswith(".gz"):
        with gzip.open(path, mode="rt", encoding="utf-8") as handle:
            consume(handle)
    else:
        with path.open(mode="r", encoding="utf-8") as handle:
            consume(handle)
    print(f"Linked {len(sequences):,} target sequence(s) from FASTA.")
    if duplicate_candidates:
        print(
            f"Resolved {duplicate_candidates:,} duplicate FASTA candidate(s) by "
            "sequence length, then sequencing-platform priority."
        )
    return sequences, annotations


def clean_dna_sequence(value: object) -> Optional[str]:
    if pd.isna(value):
        return None
    sequence = re.sub(r"\s+", "", str(value)).upper()
    return sequence or None


def write_fasta(records, path: Path, line_width: int = 80) -> None:
    """Write mapping items or an iterable of (header, sequence) pairs."""
    items = records.items() if hasattr(records, "items") else records
    with path.open("w", encoding="utf-8") as handle:
        for sequence_id, sequence in items:
            handle.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), line_width):
                handle.write(sequence[start : start + line_width] + "\n")
    print(f"Wrote {path}")


def sequence_outputs(
    df: pd.DataFrame, fasta_path: Optional[Path], output_dir: Path,
    output_prefix: str = "helotiales",
    fasta_data: Optional[tuple[dict[str, str], dict[str, dict[str, str]]]] = None,
) -> pd.DataFrame:
    """Link input/table sequences, export unique FASTA, and summarize variants."""
    target_ids = set(df["sequence_id"].dropna().astype(str))
    if fasta_data is not None:
        fasta_sequences, fasta_annotations = fasta_data
    elif fasta_path is not None:
        fasta_sequences, fasta_annotations = read_fasta(fasta_path, target_ids)
    else:
        fasta_sequences = {}
        fasta_annotations = {}
    selected: dict[str, str] = {}

    variants = df[df["sequence_id"].notna()].copy()
    for sequence_id, group in variants.groupby("sequence_id", sort=False):
        sequence_id = str(sequence_id)
        candidate = select_sh_sequence_candidate(group, fasta_sequences, fasta_annotations)
        if candidate is not None:
            selected[sequence_id] = str(candidate["sequence"])

    def first_value(series: pd.Series):
        non_null = series.dropna()
        return non_null.iloc[0] if not non_null.empty else pd.NA

    if variants.empty:
        summary = pd.DataFrame(
            columns=[
                "sequence_id", "taxonomy", "genus", "sh_id", "n_samples",
                "n_countries", "sequence_length",
            ]
        )
    else:
        summary = (
            variants.groupby("sequence_id", sort=False)
            .agg(
                taxonomy=("taxonomy", first_value),
                genus=("genus", first_value),
                sh_id=("sh_id", first_value),
                n_samples=("sample_id", safe_nunique),
                n_countries=("country", safe_nunique),
            )
            .reset_index()
        )
        summary["sequence_length"] = summary["sequence_id"].astype(str).map(
            lambda value: len(selected[value]) if value in selected else pd.NA
        )

    if selected:
        write_fasta(selected, output_dir / f"{output_prefix}_sequences.fasta")
    else:
        print("No linkable DNA sequences found; no FASTA output was written.")
    return summary


def fasta_safe(value: object) -> str:
    text = clean_metadata_text(value) or "Unknown"
    text = re.sub(r"\s+", "_", text)
    return re.sub(r"[|>\r\n]+", "_", text)


def first_useful(series: pd.Series) -> object:
    for value in series:
        cleaned = clean_metadata_text(value)
        if cleaned:
            return cleaned
    return pd.NA


def normalized_its_region(value: object) -> str:
    text = normalized_name(clean_metadata_text(value) or "")
    if text == "its1":
        return "ITS1"
    if text == "its2":
        return "ITS2"
    if text in {"itsboth", "its1its2", "both"}:
        return "ITSboth"
    return "Unknown"


def select_sh_sequence_candidate(
    sh_group: pd.DataFrame,
    sequences: dict[str, str],
    fasta_annotations: dict[str, dict[str, str]],
) -> Optional[dict[str, object]]:
    """Select one sequence for an SH by length, platform, then stable identifier."""
    candidates: list[dict[str, object]] = []
    sequence_rows = sh_group[sh_group["sequence_id"].notna()].copy()
    for sequence_id, sequence_group in sequence_rows.groupby("sequence_id", sort=False):
        sequence_id = str(sequence_id)
        fasta_sequence = sequences.get(sequence_id)
        annotation = fasta_annotations.get(sequence_id, {})
        if fasta_sequence:
            metadata_platforms = list(sequence_group["sequencing_platform"])
            header_platform = annotation.get("sequencing_platform", "Unknown")
            platform_name = best_sequencing_platform([*metadata_platforms, header_platform])
            metadata_available = any(clean_metadata_text(value) for value in metadata_platforms)
            platform_source = (
                "sample_metadata"
                if metadata_available
                else "FASTA_header"
                if normalized_sequencing_platform(header_platform) != "Unknown"
                else "unavailable"
            )
            candidates.append({
                "sequence_id": sequence_id,
                "sequence": fasta_sequence,
                "platform": platform_name,
                "platform_source": platform_source,
                "source": "FASTA",
                "accession": annotation.get("accession", sequence_id),
                "annotation": annotation,
                "represented_candidates": int(annotation.get("fasta_candidate_count", "1")),
            })

        table_work = sequence_group.copy()
        table_work["_candidate_sequence"] = table_work["sequence"].map(clean_dna_sequence)
        table_work = table_work[table_work["_candidate_sequence"].notna()]
        for table_sequence, rows in table_work.groupby("_candidate_sequence", sort=False):
            platform_name = best_sequencing_platform(rows["sequencing_platform"])
            candidates.append({
                "sequence_id": sequence_id,
                "sequence": str(table_sequence),
                "platform": platform_name,
                "platform_source": "sample_metadata",
                "source": "occurrence_table",
                "accession": sequence_id,
                "annotation": annotation,
                "represented_candidates": 1,
            })

    if not candidates:
        return None
    return min(
        candidates,
        key=lambda candidate: sequence_selection_key(
            str(candidate["sequence"]),
            candidate["platform"],
            f"{candidate['sequence_id']}|{candidate['accession']}|{candidate['source']}",
        ),
    ) | {"candidate_sequence_count": sum(
        int(candidate["represented_candidates"]) for candidate in candidates
    )}


def root_host_fasta_outputs(
    root: pd.DataFrame,
    confidence_occurrence_subsets: dict[str, pd.DataFrame],
    fasta_data: tuple[dict[str, str], dict[str, dict[str, str]]],
    output_dir: Path,
    output_prefix: str,
) -> dict[str, pd.DataFrame]:
    """Write confidence x ITS FASTAs with identical representatives per ITS set."""
    sequences, fasta_annotations = fasta_data
    work = root.copy()
    required_confidence_labels = ("medium_or_high",)
    missing_confidence = [
        label for label in required_confidence_labels
        if label not in confidence_occurrence_subsets
    ]
    if missing_confidence:
        raise ValueError(
            "Missing confidence occurrence subset(s) for FASTA output: "
            + ", ".join(missing_confidence)
        )
    confidence_count_maps = {
        label: (
            confidence_occurrence_subsets[label]
            .dropna(subset=["sh_id"])
            .groupby("sh_id", sort=False)
            .size()
            .astype(int)
            .to_dict()
        )
        for label in required_confidence_labels
    }
    work["_normalized_its_region"] = work["its_region"].map(normalized_its_region)
    summaries: dict[str, pd.DataFrame] = {}
    summary_columns = [
        "sequence_id", "sh_id", "species", "genus",
        "confidence_subset", "its_subset", "host_candidates", "countries",
        "its_regions", "sample_ids",
        "n_sh_sample_ids", "n_unique_occurrences", "n_occurrence_records",
        "n_member_sequence_ids", "sample_count_basis",
        "unique_occurrence_count_basis",
        "sequence_length", "selected_sequencing_platform", "sequencing_platform_source",
        "selected_sequence_source", "selected_accession", "candidate_sequence_count",
        "selection_rule", "fasta_header",
    ]
    region_subsets = {
        "all": work,
        "its1_only": work[work["_normalized_its_region"] == "ITS1"],
        "its2_only": work[work["_normalized_its_region"] == "ITS2"],
        "its1_and_its2": work[work["_normalized_its_region"] == "ITSboth"],
    }
    for region_label, region_subset in region_subsets.items():
        subset = region_subset.copy()
        subset["_output_sh_id"] = subset["sh_id"].map(clean_metadata_text)
        subset["_output_sh_id"] = subset["_output_sh_id"].fillna(
            subset["sequence_id"].map(clean_metadata_text)
        )
        representative_rows: list[dict[str, object]] = []
        for sh_id, group in subset[subset["_output_sh_id"].notna()].groupby(
            "_output_sh_id", sort=False
        ):
            sh_id = str(sh_id)
            selected = select_sh_sequence_candidate(group, sequences, fasta_annotations)
            if selected is None:
                continue
            fasta_annotation = selected["annotation"]
            species = (
                clean_metadata_text(first_useful(group["species"]))
                or fasta_annotation.get("species", "Unknown")
            )
            genus = (
                clean_metadata_text(first_useful(group["genus"]))
                or fasta_annotation.get("genus", "Unknown")
            )
            representative_rows.append({
                "sequence_id": str(selected["sequence_id"]),
                "sh_id": sh_id,
                "species": species,
                "genus": genus,
                "sequence": str(selected["sequence"]),
                "host_candidates": join_unique_text(group["host_candidate"]),
                "countries": join_unique_text(group["country"]),
                "its_regions": join_unique_text(group["its_region"]),
                "sample_ids": join_unique_text(group["sample_id"]),
                "n_sh_sample_ids": int(group["sample_id"].nunique(dropna=True)),
                "n_occurrence_records": len(group),
                "n_member_sequence_ids": int(group["sequence_id"].nunique(dropna=True)),
                "selected": selected,
            })

        for confidence_label in required_confidence_labels:
            output_label = region_label
            records: list[tuple[str, str]] = []
            rows: list[dict[str, object]] = []
            count_map = confidence_count_maps[confidence_label]
            for representative in representative_rows:
                sh_id = str(representative["sh_id"])
                species = representative["species"]
                sequence = str(representative["sequence"])
                selected = representative["selected"]
                n_unique_occurrences = int(count_map.get(sh_id, 0))
                header = "|".join(
                    fasta_safe(value)
                    for value in [
                        sh_id,
                        species,
                        f"n_unique_occurrences={n_unique_occurrences}",
                    ]
                )
                records.append((header, sequence))
                rows.append({
                    "sequence_id": representative["sequence_id"],
                    "sh_id": sh_id,
                    "species": species,
                    "genus": representative["genus"],
                    "confidence_subset": confidence_label,
                    "its_subset": region_label,
                    "host_candidates": representative["host_candidates"],
                    "countries": representative["countries"],
                    "its_regions": representative["its_regions"],
                    "sample_ids": representative["sample_ids"],
                    "n_sh_sample_ids": representative["n_sh_sample_ids"],
                    "n_unique_occurrences": n_unique_occurrences,
                    "n_occurrence_records": representative["n_occurrence_records"],
                    "n_member_sequence_ids": representative["n_member_sequence_ids"],
                    "sample_count_basis": (
                        "unique_sample_id_across_all_member_sequences_in_SH"
                    ),
                    "unique_occurrence_count_basis": (
                        "rows_for_SH_in_root_occurrences_with_hosts"
                    ),
                    "sequence_length": sequence_length(sequence),
                    "selected_sequencing_platform": selected["platform"],
                    "sequencing_platform_source": selected["platform_source"],
                    "selected_sequence_source": selected["source"],
                    "selected_accession": selected["accession"],
                    "candidate_sequence_count": selected["candidate_sequence_count"],
                    "selection_rule": (
                        "longest_then_platform:Sanger>PacBio>Illumina>DNBSEQ>Ion Torrent>"
                        "454 Roche>Other>Unknown"
                    ),
                    "fasta_header": header,
                })
            fasta_path = (
                output_dir / f"{output_prefix}_root_host_sequences_{output_label}.fasta"
            )
            write_fasta(records, fasta_path)
            if not records:
                print(f"FASTA subset is empty: {output_label}")
            summaries[output_label] = pd.DataFrame(rows, columns=summary_columns)
    return summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--occurrences", type=Path, default=OCCURRENCE_TABLE)
    parser.add_argument("--metadata", type=Path, default=SAMPLE_METADATA_TABLE)
    parser.add_argument(
        "--taxonomy", type=Path, default=TAXONOMY_TABLE,
        help="SH-to-taxonomy table; required for a wide SH abundance matrix unless --sh-list is used.",
    )
    parser.add_argument(
        "--sh-list", type=Path, default=SH_LIST_FILE,
        help="Optional text/table containing explicit SH IDs to extract.",
    )
    parser.add_argument(
        "--target-taxon", action="append", dest="target_taxa",
        help="Taxon name to select; repeat for multiple taxa (default: Helotiales).",
    )
    parser.add_argument("--fasta", type=Path, default=FASTA_FILE)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument(
        "--figure-font-scale",
        type=float,
        default=1.5,
        help="Multiplier for every figure text size (default: 1.5).",
    )
    parser.add_argument(
        "--world-geojson",
        type=Path,
        default=WORLD_BOUNDARIES_FILE,
        help=(
            "GeoJSON country/land boundaries for the sampling map "
            "(default: bundled Natural Earth asset)."
        ),
    )
    parser.add_argument("--minimum-reads", type=float, default=MINIMUM_READ_COUNT)
    parser.add_argument(
        "--matrix-chunk-size", type=int, default=WIDE_MATRIX_CHUNK_SIZE,
        help="Number of sample rows per chunk for a wide abundance matrix.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=0,
        help=(
            "Reserved parallel-worker setting. Use 0 for all available logical CPUs, "
            "or a positive integer to set it explicitly."
        ),
    )
    parser.add_argument(
        "--read-count-weighting",
        action="store_true",
        help="Set analysis_weight to reads instead of sample presence (not recommended globally).",
    )
    args = parser.parse_args()
    if not math.isfinite(args.figure_font_scale) or args.figure_font_scale <= 0:
        parser.error("--figure-font-scale must be a finite number greater than 0")
    return args


def clean_target_taxa(values: Optional[list[str]]) -> list[str]:
    """Accept repeated options and comma-separated names while preserving order."""
    source = values if values else TARGET_TAXA
    targets = []
    for value in source:
        targets.extend(part.strip() for part in str(value).split(",") if part.strip())
    return list(dict.fromkeys(targets))


def output_slug(targets: list[str]) -> str:
    if len(targets) == 1:
        slug = re.sub(r"[^a-z0-9]+", "_", targets[0].casefold()).strip("_")
        return slug or "target_taxon"
    return "selected_taxa"


def discover_taxonomy_path(
    supplied: Optional[Path], occurrence_path: Path
) -> Optional[Path]:
    """Auto-resolve the documented placeholder to one local SH taxonomy file."""
    if supplied is not None and supplied.exists():
        return supplied
    parent = occurrence_path.resolve().parent
    candidates = sorted(
        set(parent.glob("*SHs.tax*"))
        | set(parent.glob("*SH_taxonomy*.tsv"))
    )
    is_placeholder = supplied is not None and supplied.name == "SH_taxonomy.tsv"
    if len(candidates) == 1 and (supplied is None or is_placeholder):
        print(f"Auto-detected SH taxonomy source: {candidates[0]}")
        return candidates[0]
    return supplied


def discover_fasta_path(
    supplied: Optional[Path],
    occurrence_path: Path,
    taxonomy_path: Optional[Path] = None,
) -> Optional[Path]:
    """Use an explicit FASTA path or auto-detect one UNITE archive beside inputs."""
    if supplied is not None:
        return supplied
    parent = occurrence_path.resolve().parent
    candidates = sorted(
        set(parent.glob("sh_general_release*.tgz"))
        | set(parent.glob("sh_general_release*.tar.gz"))
        | set(parent.glob("sh_general_release*.fasta"))
    )
    if len(candidates) == 1:
        print(f"Auto-detected UNITE FASTA source: {candidates[0]}")
        return candidates[0]
    if len(candidates) > 1:
        taxonomy_date = None
        if taxonomy_path is not None:
            match = re.search(r"\d{2}\.\d{2}\.\d{4}", taxonomy_path.name)
            taxonomy_date = match.group(0) if match else None
        matching = [path for path in candidates if taxonomy_date and taxonomy_date in path.name]
        if len(matching) == 1:
            print(
                f"Auto-selected UNITE FASTA matching taxonomy release "
                f"{taxonomy_date}: {matching[0]}"
            )
            return matching[0]
        warnings.warn("Multiple UNITE FASTA sources found; specify the intended file with --fasta.")
    return None


def main() -> int:
    args = parse_args()
    set_figure_font_scale(args.figure_font_scale)
    effective_jobs = resolve_parallel_jobs(args.jobs)
    targets = clean_target_taxa(args.target_taxa)
    target_label = ", ".join(targets)
    prefix = output_slug(targets)
    taxonomy_path = discover_taxonomy_path(args.taxonomy, args.occurrences)
    fasta_path = discover_fasta_path(args.fasta, args.occurrences, taxonomy_path)
    if args.matrix_chunk_size < 1:
        raise ValueError("--matrix-chunk-size must be at least 1.")
    output_dir = args.output.resolve()
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    initialize_run_log(output_dir, args, effective_jobs)
    version_path = output_dir / "runtime_versions.tsv"
    save_tsv(
        pd.DataFrame(
            runtime_version_records(Path(__file__).resolve(), args.jobs, effective_jobs)
        ),
        version_path,
    )
    print(
        f"Parallel workers: requested={args.jobs}, resolved={effective_jobs} "
        "(recorded for reproducibility; no current stage needs worker processes)."
    )

    layout, occurrence_columns, occurrence_delimiter = detect_occurrence_layout(args.occurrences)
    print(
        f"Occurrence input layout: {layout} "
        f"({len(occurrence_columns):,} columns; targets: {target_label})"
    )
    metadata_raw = load_table(args.metadata)
    metadata_mapping = detect_columns(metadata_raw)
    inspect_table("Sample metadata table", metadata_raw, metadata_mapping)
    if metadata_mapping["sample_id"] is None:
        raise ValueError(
            "No sample-ID column was detected in the metadata table. Add its name "
            "to COLUMN_ALIASES['sample_id']."
        )
    metadata = standardize_table(metadata_raw, metadata_mapping)

    if layout == "wide_sh_matrix":
        if args.sh_list is not None:
            selection = load_explicit_sh_list(args.sh_list)
            selection = annotate_explicit_shs(selection, taxonomy_path)
            match_summary = pd.DataFrame([{
                "criterion": "explicit_sh_list",
                "n_records": len(selection),
            }])
        elif taxonomy_path is not None:
            if not taxonomy_path.exists():
                raise FileNotFoundError(
                    f"The taxonomy file supplied with --taxonomy does not exist: "
                    f"{taxonomy_path}. 'data/SH_taxonomy.tsv' in the usage example is a "
                    "placeholder, not a file included with GlobalFungi. Download/provide "
                    "a real SH taxonomy table, or use --sh-list instead."
                )
            selection, match_summary = select_shs_from_taxonomy(taxonomy_path, targets)
        else:
            raise ValueError(
                "The occurrence input is a wide sample-by-SH matrix. Supply an SH taxonomy "
                "mapping with --taxonomy, or explicit target SH IDs with --sh-list. The "
                "abundance matrix alone contains no taxonomic assignments."
            )
        resolved = resolve_matrix_sh_columns(occurrence_columns, selection)
        selected_occurrences = extract_wide_sh_occurrences(
            args.occurrences,
            occurrence_columns,
            occurrence_delimiter,
            resolved,
            args.minimum_reads,
            args.matrix_chunk_size,
        )
        occurrence_mapping = detect_columns(selected_occurrences)
        selected_occurrences = standardize_table(selected_occurrences, occurrence_mapping)
    else:
        occurrence_raw = load_table(args.occurrences)
        occurrence_mapping = detect_columns(occurrence_raw)
        inspect_table("Occurrence table", occurrence_raw, occurrence_mapping)
        validate_required_keys(occurrence_mapping, metadata_mapping)
        occurrence = standardize_table(occurrence_raw, occurrence_mapping)
        occurrence = filter_minimum_reads(
            occurrence, occurrence_mapping["reads"], args.minimum_reads
        )
        selected, reasons, match_summary = match_target_taxa(occurrence, targets)
        selected_occurrences = occurrence.loc[selected].copy()
        selected_occurrences["target_taxon_match_reason"] = reasons.loc[selected]
        selected_occurrences["helotiales_match_reason"] = reasons.loc[selected]

    print("\nTarget-taxon selection criteria (criteria can overlap):")
    print(match_summary.to_string(index=False))
    save_tsv(match_summary, output_dir / f"{prefix}_match_summary.tsv")
    if selected_occurrences.empty:
        print(
            "\nNo records matched. Check the taxonomy mapping, target spelling, and aliases.",
            file=sys.stderr,
        )
        return 2

    merged = merge_occurrences_metadata(selected_occurrences, metadata)
    merged = apply_host_metadata_priority(merged)
    merged["presence"] = 1
    use_presence = USE_PRESENCE_ABSENCE and not args.read_count_weighting
    if use_presence:
        merged["analysis_weight"] = 1.0
        merged["analysis_weight_type"] = "presence_absence"
    else:
        merged["analysis_weight"] = merged["reads"]
        merged["analysis_weight_type"] = "read_count"
        warnings.warn(
            "Read-count weighting is enabled. Cross-study read counts are not directly comparable."
        )

    dataset_stats = report_merged_data(merged, target_label)
    save_tsv(merged, output_dir / f"{prefix}_occurrences_merged.tsv")
    save_tsv(dataset_stats, output_dir / "dataset_summary.tsv")

    save_tsv(
        geographic_ranges_by_genus(merged),
        output_dir / "geographic_ranges_by_genus.tsv",
    )

    root_records, root_occurrences = build_root_host_tables(merged)
    occurrence_confidence_subsets: dict[str, pd.DataFrame] = {}
    print(
        f"\nRoot-only subset: {len(root_occurrences):,} unique occurrences, "
        f"{root_occurrences['sample_id'].nunique(dropna=True) if not root_occurrences.empty else 0:,} samples"
    )
    if not root_occurrences.empty:
        root_occurrences = validate_root_occurrences(root_occurrences)
        occurrence_confidence = (
            root_occurrences["host_confidence"].astype("string").str.casefold()
        )
        occurrence_confidence_subsets = {
            "medium_or_high": root_occurrences[
                occurrence_confidence.isin(["medium", "high"])
            ].copy(),
        }
        root_occurrences = occurrence_confidence_subsets["medium_or_high"].copy()
        root_occurrence_path = output_dir / f"{prefix}_root_occurrences_with_hosts.tsv"
        save_tsv(root_occurrences, root_occurrence_path)
        verify_tsv_data_lines_start_with_sh(root_occurrence_path)
        unique_occurrence_validation = pd.DataFrame([
            {
                "confidence_subset": label,
                "n_unique_occurrences": len(subset),
                "n_source_sample_sh_occurrences": int(
                    subset["n_sample_sh_occurrences"].sum()
                ),
                "n_sh_ids": int(subset["sh_id"].nunique(dropna=True)),
                "n_representative_sample_ids": int(
                    subset["sample_id"].nunique(dropna=True)
                ),
            }
            for label, subset in occurrence_confidence_subsets.items()
        ])
        save_tsv(
            unique_occurrence_validation,
            output_dir / f"{prefix}_unique_occurrence_count_validation.tsv",
        )
    fasta_data = None
    if fasta_path is not None:
        target_sequence_ids = set(merged["sequence_id"].dropna().astype(str))
        fasta_data = read_fasta(fasta_path, target_sequence_ids)
        linked_ids = set(fasta_data[0])
        linkage = pd.DataFrame({"sequence_id": sorted(target_sequence_ids)})
        linkage["fasta_sequence_found"] = linkage["sequence_id"].isin(linked_ids)
        save_tsv(linkage, output_dir / f"{prefix}_fasta_linkage_report.tsv")
        if not linked_ids and target_sequence_ids:
            warnings.warn(
                "No selected SH IDs were found in the supplied FASTA source. This usually "
                "means that the UNITE FASTA release does not match the SH release used by "
                "GlobalFungi. Use the corresponding UNITE release or GlobalFungi-exported sequences."
            )

    sequence_summary = sequence_outputs(
        merged, fasta_path, output_dir, prefix, fasta_data=fasta_data
    )
    save_tsv(sequence_summary, output_dir / f"{prefix}_sequence_summary.tsv")
    if not root_records.empty and fasta_data is not None:
        root_sequence_summaries = root_host_fasta_outputs(
            root_records,
            occurrence_confidence_subsets,
            fasta_data,
            output_dir,
            prefix,
        )
        for label, root_sequence_summary in root_sequence_summaries.items():
            save_tsv(
                root_sequence_summary,
                output_dir / f"{prefix}_root_host_sequence_summary_{label}.tsv",
            )

    plot_world_points(
        merged,
        figure_dir / "world_map_occurrences.pdf",
        args.world_geojson,
        f"{target_label}: all Helotiales samples",
        occurrence_label="Occurrences",
    )
    if not root_occurrences.empty:
        plot_world_points(
            root_occurrences,
            figure_dir / "world_map_root_unique_occurrences.pdf",
            args.world_geojson,
            f"{target_label}: medium-or-higher root unique occurrences",
            occurrence_label="Unique occurrences",
        )
    plot_latitude_histogram(merged, figure_dir / "latitude_histogram.pdf")

    print(
        "\nFinished. Presence/absence is the default analysis weight. Runtime and package "
        f"versions are recorded in {version_path.name}."
    )
    print(f"Output directory: {output_dir}")
    return 0


if __name__ == "__main__":
    exit_code = 1
    try:
        exit_code = main()
    except BaseException:
        finalize_run_log("failed", exit_code)
        raise
    else:
        finalize_run_log("completed" if exit_code == 0 else "nonzero_exit", exit_code)
        raise SystemExit(exit_code)
