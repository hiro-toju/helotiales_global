#!/usr/bin/env python3
"""Draw occurrence-record maps for abundant Helotiales SHs.

This script is a small, standalone extraction of the SH spatial-map output that
used to be embedded in the 2_Specificity workflow.  It reads the unique
occurrence table from 1_GlobalFungi, ranks SHs by unique sampling-unit counts,
and writes one PDF map per selected SH.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import shlex
import sys
import time
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
UNKNOWN_LABELS = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "unknown",
    "unidentified",
    "unclassified",
    "unresolved",
}
INVALID_PLANT_GENERA = {
    "viridiplantae",
    "plantae",
    "tracheophyta",
    "streptophyta",
    "embryophyta",
    "vascular plants",
    "plant",
    "plants",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create PDF occurrence maps for abundant SHs from the "
            "1_GlobalFungi unique-occurrence table."
        )
    )
    parser.add_argument(
        "--globalfungi-results",
        type=Path,
        default=PROJECT_DIR / "1_GlobalFungi" / "results_helotiales",
        help=(
            "Directory containing helotiales_root_occurrences_with_hosts.tsv "
            "(default: ../1_GlobalFungi/results_helotiales)."
        ),
    )
    parser.add_argument(
        "--occurrence-file",
        type=Path,
        default=None,
        help=(
            "Optional direct path to a unique-occurrence TSV. If omitted, "
            "--globalfungi-results/helotiales_root_occurrences_with_hosts.tsv is used."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Output root directory (default: outputs).",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help=(
            "Output subfolder name. If omitted, it is generated from the "
            "selection settings, for example sh50_maps."
        ),
    )
    parser.add_argument(
        "--sh-selection",
        choices=("rank", "min-count", "list"),
        default="rank",
        help=(
            "How to select SHs: top rank, minimum unique-occurrence count, "
            "or explicit --selected-sh-file (default: rank)."
        ),
    )
    parser.add_argument(
        "--sh-rank-max",
        type=int,
        default=50,
        help="Number of top-ranked SHs to map when --sh-selection rank (default: 50).",
    )
    parser.add_argument(
        "--sh-min-count",
        type=int,
        default=1,
        help=(
            "Minimum unique-occurrence count when --sh-selection min-count "
            "(default: 1)."
        ),
    )
    parser.add_argument(
        "--selected-sh-file",
        type=Path,
        default=None,
        help="Optional TSV/CSV/plain-text file containing SH IDs to map.",
    )
    parser.add_argument(
        "--selected-sh-column",
        default="sh_id",
        help="Column containing SH IDs in --selected-sh-file (default: sh_id).",
    )
    parser.add_argument(
        "--plant-taxonomy",
        type=Path,
        default=SCRIPT_DIR / "data" / "plant_taxonomy.tsv",
        help=(
            "Cached plant taxonomy TSV used to remove non-genus plant labels "
            "(default: data/plant_taxonomy.tsv)."
        ),
    )
    parser.add_argument(
        "--keep-uncertain-hosts",
        action="store_true",
        help=(
            "Keep occurrences whose host label is ambiguous or not resolved as a "
            "plant genus. By default these are excluded to mimic 2_Specificity maps."
        ),
    )
    parser.add_argument(
        "--world-geojson",
        type=Path,
        default=SCRIPT_DIR / "assets" / "ne_110m_admin_0_countries.geojson",
        help="Natural Earth country-boundary GeoJSON for land outlines.",
    )
    parser.add_argument(
        "--figure-font-scale",
        type=float,
        default=1.5,
        help="Multiplier for figure font sizes (default: 1.5).",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=13.0,
        help="Map point size (default: 13).",
    )
    parser.add_argument(
        "--point-alpha",
        type=float,
        default=0.65,
        help="Map point transparency (default: 0.65).",
    )
    return parser.parse_args()


def portable_path(path: Path, base: Path) -> str:
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except ValueError:
        try:
            return str(path.resolve().relative_to(Path.cwd().resolve()))
        except ValueError:
            return str(path)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "unnamed"


def is_unknown(value: object) -> bool:
    text = "" if pd.isna(value) else str(value).strip()
    return text.casefold() in UNKNOWN_LABELS


def host_plant_genus(value: object) -> str:
    """Return a conservative plant-genus candidate from a host label."""
    if pd.isna(value):
        return "Unknown"
    text = str(value).strip()
    if is_unknown(text):
        return "Unknown"
    if re.search(r"[,;/|]", text) or re.search(
        r"\b(?:and|mixed|multiple|various|several)\b", text, flags=re.I
    ):
        return "Unknown"
    first = text.split()[0].strip("()[]{}.,;:'\"")
    if not re.fullmatch(r"[A-Z][A-Za-z-]+", first):
        return "Unknown"
    if first.casefold().endswith(("aceae", "eae", "oideae", "phyta", "opsida")):
        return "Unknown"
    return first


def load_plant_taxonomy(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")


def plant_genus_validity(plant_genera: pd.Series, taxonomy: pd.DataFrame | None) -> pd.Series:
    valid = ~plant_genera.map(is_unknown)
    valid &= ~plant_genera.astype(str).str.casefold().isin(INVALID_PLANT_GENERA)
    if taxonomy is None:
        return valid
    tax = taxonomy.set_index("input_name", drop=False)
    checked = []
    for genus in plant_genera.astype(str):
        if genus not in tax.index:
            checked.append(False)
            continue
        record = tax.loc[genus]
        accepted_rank = str(record.get("accepted_rank", "")).strip().upper()
        resolved_genus = str(record.get("genus", "")).strip()
        checked.append(
            accepted_rank == "GENUS"
            and bool(resolved_genus)
            and resolved_genus.casefold() not in INVALID_PLANT_GENERA
        )
    return valid & pd.Series(checked, index=plant_genera.index)


def read_occurrences(path: Path, plant_taxonomy: Path, keep_uncertain_hosts: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"sh_id", "latitude", "longitude", "host_candidate"}
    header = pd.read_csv(path, sep="\t", nrows=0, encoding="utf-8-sig").columns
    missing = sorted(required - set(header))
    if missing:
        raise ValueError(f"Occurrence file is missing required columns: {', '.join(missing)}")
    data = pd.read_csv(path, sep="\t", dtype=str, encoding="utf-8-sig", low_memory=False)
    data["sh_id"] = data["sh_id"].astype("string").str.strip()
    data["latitude"] = pd.to_numeric(data["latitude"], errors="coerce")
    data["longitude"] = pd.to_numeric(data["longitude"], errors="coerce")
    data["occurrence_count"] = pd.to_numeric(
        data.get("occurrence_count", 1), errors="coerce"
    ).fillna(1).astype(int)
    data["plant_genus_candidate"] = data["host_candidate"].map(host_plant_genus)
    data["map_record_id"] = np.arange(1, len(data) + 1)

    reasons = pd.Series("", index=data.index, dtype="string")

    def add_reason(mask: pd.Series, reason: str) -> None:
        current = reasons.loc[mask].fillna("")
        reasons.loc[mask] = np.where(current.eq(""), reason, current + ";" + reason)

    add_reason(data["sh_id"].isna() | ~data["sh_id"].str.fullmatch(r"SH\d+(?:\.\d+FU)?", na=False), "invalid_sh_id")
    add_reason(data["latitude"].isna() | data["longitude"].isna(), "missing_coordinates")
    add_reason(~data["latitude"].between(-90, 90) | ~data["longitude"].between(-180, 180), "coordinates_out_of_range")

    if not keep_uncertain_hosts:
        taxonomy = load_plant_taxonomy(plant_taxonomy)
        valid_host = plant_genus_validity(data["plant_genus_candidate"], taxonomy)
        add_reason(~valid_host, "uncertain_or_non_genus_host_label")

    excluded = data.loc[reasons.ne("")].copy()
    excluded.insert(0, "exclusion_reason", reasons.loc[excluded.index])
    kept = data.loc[reasons.eq("")].copy()
    sampling_cols = ["latitude", "longitude", "host_candidate"]
    kept["sampling_unit_key"] = (
        kept[sampling_cols].astype(str).agg("|".join, axis=1)
    )
    return kept.reset_index(drop=True), excluded.reset_index(drop=True)


def rank_shs(occurrences: pd.DataFrame) -> pd.DataFrame:
    unit_table = occurrences.drop_duplicates(["sh_id", "sampling_unit_key"])
    grouped = (
        unit_table.groupby("sh_id", sort=False)
        .agg(
            n_unique_occurrences=("sampling_unit_key", "nunique"),
            n_sampling_sites=("latitude", lambda x: len(
                unit_table.loc[x.index, ["latitude", "longitude"]].drop_duplicates()
            )),
            n_sample_ids=("sample_id", "nunique") if "sample_id" in unit_table else ("sampling_unit_key", "size"),
        )
        .reset_index()
    )
    tax_cols = [c for c in ["species", "genus", "family", "order", "class"] if c in occurrences.columns]
    if tax_cols:
        taxonomy = (
            occurrences[["sh_id", *tax_cols]]
            .replace({"": pd.NA, "nan": pd.NA})
            .dropna(how="all", subset=tax_cols)
            .drop_duplicates("sh_id")
        )
        grouped = grouped.merge(taxonomy, on="sh_id", how="left")
    grouped = grouped.sort_values(
        ["n_unique_occurrences", "sh_id"], ascending=[False, True]
    ).reset_index(drop=True)
    grouped.insert(0, "rank", np.arange(1, len(grouped) + 1))
    return grouped


def read_selected_shs(path: Path, column: str) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Selected SH file does not exist: {path}")
    if path.suffix.lower() in {".tsv", ".txt", ".csv"}:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        table = pd.read_csv(path, sep=sep, dtype=str, comment="#")
        if column in table.columns:
            return table[column].dropna().astype(str).str.strip().loc[lambda s: s.ne("")].tolist()
        if table.shape[1] == 1:
            return table.iloc[:, 0].dropna().astype(str).str.strip().loc[lambda s: s.ne("")].tolist()
        raise ValueError(f"Column {column!r} was not found in {path}")
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def select_shs(ranking: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    ranking = ranking.copy()
    ranking["selected_for_map"] = "not_selected"
    if args.sh_selection == "list":
        if args.selected_sh_file is None:
            raise ValueError("--sh-selection list requires --selected-sh-file")
        wanted = list(dict.fromkeys(read_selected_shs(args.selected_sh_file, args.selected_sh_column)))
        order = {sh: i for i, sh in enumerate(wanted)}
        selected = ranking.loc[ranking["sh_id"].isin(order)].copy()
        selected["_list_order"] = selected["sh_id"].map(order)
        selected = selected.sort_values("_list_order").drop(columns="_list_order")
    elif args.sh_selection == "min-count":
        selected = ranking.loc[ranking["n_unique_occurrences"].ge(args.sh_min_count)].copy()
    else:
        selected = ranking.head(args.sh_rank_max).copy()
    ranking.loc[ranking["sh_id"].isin(selected["sh_id"]), "selected_for_map"] = "selected"
    return selected, ranking


def draw_world_boundaries(ax: plt.Axes, geojson_path: Path) -> None:
    if not geojson_path.exists():
        raise FileNotFoundError(f"World-boundary GeoJSON does not exist: {geojson_path}")
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


def apply_font_scale(fig: plt.Figure, scale: float) -> None:
    for obj in list(fig.findobj(match=plt.Text)):
        obj.set_fontsize(obj.get_fontsize() * scale)


def map_one_sh(
    occurrences: pd.DataFrame,
    sh_id: str,
    rank: int,
    output_dir: Path,
    geojson_path: Path,
    font_scale: float,
    point_size: float,
    point_alpha: float,
) -> Path:
    subset = occurrences.loc[occurrences["sh_id"].eq(sh_id)].copy()
    n_units = int(subset.drop_duplicates(["sh_id", "sampling_unit_key"]).shape[0])
    points = subset[["latitude", "longitude"]].drop_duplicates()
    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.set_facecolor("#f7fbfc")
    draw_world_boundaries(ax, geojson_path)
    if not points.empty:
        ax.scatter(
            points["longitude"],
            points["latitude"],
            s=point_size,
            alpha=point_alpha,
            color="#8B2F67",
            edgecolors="none",
            zorder=2,
        )
    ax.set(xlim=(-180, 180), ylim=(-90, 90), xlabel="Longitude", ylabel="Latitude")
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.grid(alpha=0.18, lw=0.4)
    ax.set_title(
        f"{sh_id}: unique occurrence records "
        f"(n = {n_units:,}; mapped sites = {len(points):,})"
    )
    apply_font_scale(fig, font_scale)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"{rank:03d}_{safe_name(sh_id)}_spatial_distribution.pdf"
    fig.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return out


def write_versions(path: Path) -> None:
    versions = [
        ("python", sys.version.replace("\n", " ")),
        ("python_executable", sys.executable),
        ("platform", platform.platform()),
        ("pandas", pd.__version__),
        ("numpy", np.__version__),
        ("matplotlib", matplotlib.__version__),
    ]
    pd.DataFrame(versions, columns=["package", "version"]).to_csv(
        path, sep="\t", index=False
    )


def run_name(args: argparse.Namespace) -> str:
    if args.run_name:
        return safe_name(args.run_name)
    if args.sh_selection == "min-count":
        return f"sh_min{args.sh_min_count}_maps"
    if args.sh_selection == "list":
        stem = args.selected_sh_file.stem if args.selected_sh_file else "list"
        return f"{safe_name(stem)}_maps"
    return f"sh{args.sh_rank_max}_maps"


def main() -> int:
    start = time.time()
    args = parse_args()
    plt.rcParams["font.family"] = "Arial"

    occurrence_file = (
        args.occurrence_file
        if args.occurrence_file is not None
        else args.globalfungi_results / "helotiales_root_occurrences_with_hosts.tsv"
    )
    if not occurrence_file.exists():
        raise FileNotFoundError(f"Occurrence file does not exist: {occurrence_file}")
    if not args.world_geojson.exists():
        raise FileNotFoundError(f"World-boundary GeoJSON does not exist: {args.world_geojson}")

    output_root = args.output_dir / run_name(args)
    maps_dir = output_root / "selected_SH_spatial_maps"
    output_root.mkdir(parents=True, exist_ok=True)
    log_path = output_root / "analysis.log"
    command_path = output_root / "command.txt"

    command_path.write_text(shlex.join(["python3", *sys.argv]) + "\n", encoding="utf-8")
    log_lines = [
        f"started_at_epoch\t{start:.3f}",
        f"occurrence_file\t{portable_path(occurrence_file, output_root)}",
        f"selection\t{args.sh_selection}",
        f"keep_uncertain_hosts\t{args.keep_uncertain_hosts}",
    ]

    occurrences, excluded = read_occurrences(
        occurrence_file, args.plant_taxonomy, args.keep_uncertain_hosts
    )
    if occurrences.empty:
        raise ValueError("No mappable occurrence records remain after filtering.")

    ranking = rank_shs(occurrences)
    selected, ranking = select_shs(ranking, args)
    if selected.empty:
        raise ValueError("No SHs passed the requested selection setting.")

    ranking.to_csv(output_root / "sh_unique_occurrence_ranking.tsv", sep="\t", index=False)
    selected.to_csv(output_root / "selected_shs.tsv", sep="\t", index=False)
    excluded.to_csv(output_root / "excluded_occurrences.tsv", sep="\t", index=False)
    write_versions(output_root / "package_versions.tsv")

    for map_rank, sh_id in enumerate(selected["sh_id"].astype(str), start=1):
        map_one_sh(
            occurrences,
            sh_id,
            map_rank,
            maps_dir,
            args.world_geojson,
            args.figure_font_scale,
            args.point_size,
            args.point_alpha,
        )

    elapsed = time.time() - start
    summary = pd.DataFrame([
        {"item": "input_rows", "value": int(len(occurrences) + len(excluded))},
        {"item": "kept_rows", "value": int(len(occurrences))},
        {"item": "excluded_rows", "value": int(len(excluded))},
        {"item": "ranked_shs", "value": int(len(ranking))},
        {"item": "mapped_shs", "value": int(len(selected))},
        {"item": "elapsed_seconds", "value": round(elapsed, 3)},
    ])
    summary.to_csv(output_root / "run_summary.tsv", sep="\t", index=False)
    log_lines.extend([
        f"kept_rows\t{len(occurrences)}",
        f"excluded_rows\t{len(excluded)}",
        f"mapped_shs\t{len(selected)}",
        f"elapsed_seconds\t{elapsed:.3f}",
        "status\tcompleted",
    ])
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"Mapped {len(selected):,} SHs.")
    print(f"Output directory: {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
