#!/usr/bin/env python3
"""Summarize cleaned Helotiales unique-occurrence data properties.

The script removes unique occurrences and plant labels rejected during the
4_Specificity_Lat preprocessing, then writes a cleaned occurrence table plus
category-count summaries and PDF donut charts.
"""

from __future__ import annotations

import argparse
import math
import os
import platform
import re
import shlex
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

CACHE_DIR = Path(tempfile.gettempdir()) / "helotiales_specificity_matplotlib"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DEFAULT_RUN_NAME = "sh150_family30_10000_lat20_minblock100"
DEFAULT_TOP_N = 15
UNKNOWN_TOKENS = {"", "NA", "N/A", "NAN", "NONE", "NULL", "UNKNOWN", "<NA>"}

plt.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_DIR.resolve()).as_posix()
    except Exception:
        return path.as_posix()


def is_unknown(value: object) -> bool:
    if pd.isna(value):
        return True
    return str(value).strip().upper() in UNKNOWN_TOKENS


def clean_label(value: object) -> str:
    if is_unknown(value):
        return "Unknown"
    return str(value).strip()


def host_plant_genus(value: object) -> str:
    """Mirror 2/3_Specificity host-candidate-to-plant-genus parsing."""
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a cleaned Helotiales occurrence table and data-property summaries."
    )
    parser.add_argument("--specificity-run-name", default=DEFAULT_RUN_NAME)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Unique-occurrence TSV. Defaults to the 4_Specificity_Lat data copy.",
    )
    parser.add_argument(
        "--specificity-result-dir",
        type=Path,
        default=None,
        help="Directory containing excluded_unique_occurrences.tsv and excluded_plant_columns.tsv.",
    )
    parser.add_argument(
        "--plant-taxonomy",
        type=Path,
        default=None,
        help="Plant taxonomy TSV used by the specificity workflow.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to outputs/<run>.",
    )
    parser.add_argument("--plot-top-n", type=int, default=DEFAULT_TOP_N)
    parser.add_argument(
        "--top-sh-bar-n",
        type=int,
        default=50,
        help="Number of most frequent fungal SHs to show in the SH bar plot.",
    )
    parser.add_argument("--figure-font-scale", type=float, default=1.2)
    return parser.parse_args()


def resolve_paths(args: argparse.Namespace) -> dict[str, Path]:
    run_root = (
        PROJECT_DIR
        / "4_Specificity_Lat"
        / "outputs"
        / args.specificity_run_name
    )
    stem = "helotiales_root_occurrences_with_hosts"
    input_path = args.input or run_root / "data" / f"{stem}.tsv"
    result_dir = args.specificity_result_dir or run_root / stem
    taxonomy_path = args.plant_taxonomy or run_root / "data" / "plant_taxonomy.tsv"
    output_dir = args.output_dir or (
        SCRIPT_DIR / "outputs" / args.specificity_run_name
    )
    return {
        "input": input_path,
        "result_dir": result_dir,
        "taxonomy": taxonomy_path,
        "output_dir": output_dir,
        "stem": Path(stem),
    }


def require_files(paths: dict[str, Path]) -> None:
    stem = paths["stem"].name
    required = [
        paths["input"],
        paths["taxonomy"],
        paths["result_dir"] / f"{stem}_excluded_unique_occurrences.tsv",
        paths["result_dir"] / f"{stem}_excluded_plant_columns.tsv",
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Required file(s) missing: " + ", ".join(rel(path) for path in missing)
        )


def read_taxonomy(path: Path) -> pd.DataFrame:
    taxonomy = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    required = {"input_name", "family", "order"}
    missing = required - set(taxonomy.columns)
    if missing:
        raise ValueError(f"Plant taxonomy table is missing columns: {sorted(missing)}")
    taxonomy["input_name"] = taxonomy["input_name"].astype(str).str.strip()
    return taxonomy


def clean_occurrences(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    stem = paths["stem"].name
    raw = pd.read_csv(paths["input"], sep="\t", dtype=str, low_memory=False)
    excluded_unique = pd.read_csv(
        paths["result_dir"] / f"{stem}_excluded_unique_occurrences.tsv",
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    excluded_plants = pd.read_csv(
        paths["result_dir"] / f"{stem}_excluded_plant_columns.tsv",
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    taxonomy = read_taxonomy(paths["taxonomy"])
    tax_by_input = taxonomy.drop_duplicates("input_name").set_index("input_name")

    key_columns = ["sh_id", "latitude", "longitude", "host_candidate"]
    for column in key_columns:
        if column not in raw.columns or column not in excluded_unique.columns:
            raise ValueError(f"Missing key column {column!r} in input or exclusion table.")
    raw_key = raw[key_columns].fillna("").astype(str).agg("\t".join, axis=1)
    excluded_key = set(excluded_unique[key_columns].fillna("").astype(str).agg("\t".join, axis=1))

    data = raw.loc[~raw_key.isin(excluded_key)].copy()
    data["plant_genus_for_analysis"] = data["host_candidate"].map(host_plant_genus)

    invalid_plant_labels = set(
        excluded_plants.get("input_plant_name", pd.Series(dtype=str))
        .dropna()
        .astype(str)
        .str.strip()
    )
    before_invalid = len(data)
    data = data.loc[~data["plant_genus_for_analysis"].isin(invalid_plant_labels)].copy()
    data["plant_family_for_analysis"] = (
        data["plant_genus_for_analysis"].map(tax_by_input["family"]).map(clean_label)
    )
    data["plant_order_for_analysis"] = (
        data["plant_genus_for_analysis"].map(tax_by_input["order"]).map(clean_label)
    )

    removal_summary = pd.DataFrame([
        {
            "step": "input",
            "n_unique_occurrences": len(raw),
            "n_sh": raw["sh_id"].nunique(dropna=True),
        },
        {
            "step": "after_excluded_unique_occurrence_filter",
            "n_unique_occurrences": before_invalid,
            "n_sh": data["sh_id"].nunique(dropna=True)
            + len(set(raw.loc[raw_key.isin(excluded_key), "sh_id"]) - set(data["sh_id"])),
        },
        {
            "step": "final_cleaned",
            "n_unique_occurrences": len(data),
            "n_sh": data["sh_id"].nunique(dropna=True),
        },
    ])
    removal_summary.loc[
        removal_summary["step"].eq("after_excluded_unique_occurrence_filter"),
        "n_sh",
    ] = raw.loc[~raw_key.isin(excluded_key), "sh_id"].nunique(dropna=True)
    removal_summary["n_removed_from_previous_step"] = (
        removal_summary["n_unique_occurrences"].shift(1) - removal_summary["n_unique_occurrences"]
    )
    removal_summary.loc[0, "n_removed_from_previous_step"] = 0
    removal_summary["n_removed_from_previous_step"] = (
        removal_summary["n_removed_from_previous_step"].fillna(0).astype(int)
    )
    return data.reset_index(drop=True), removal_summary


def value_counts_table(df: pd.DataFrame, dimension: str, series: pd.Series) -> pd.DataFrame:
    labels = series.map(clean_label)
    counts = labels.value_counts(dropna=False)
    out = counts.rename_axis("category").reset_index(name="n_unique_occurrences")
    out.insert(0, "dimension", dimension)
    out["percent_unique_occurrences"] = out["n_unique_occurrences"] / len(df) * 100
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def build_breakdowns(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    dimensions = {
        "fungal_sh": df["sh_id"],
        "fungal_genus": df["genus"] if "genus" in df.columns else pd.Series(["Unknown"] * len(df)),
        "fungal_family": df["family"] if "family" in df.columns else pd.Series(["Unknown"] * len(df)),
        "plant_genus": df["plant_genus_for_analysis"],
        "plant_family": df["plant_family_for_analysis"],
        "plant_order": df["plant_order_for_analysis"],
        "country": df["country"] if "country" in df.columns else pd.Series(["Unknown"] * len(df)),
        "continent": df["continent"] if "continent" in df.columns else pd.Series(["Unknown"] * len(df)),
    }
    breakdown = pd.concat(
        [value_counts_table(df, dimension, values) for dimension, values in dimensions.items()],
        ignore_index=True,
    )
    summary_rows = []
    for dimension, sub in breakdown.groupby("dimension", sort=False):
        summary_rows.append({
            "dimension": dimension,
            "n_categories_including_unknown": sub["category"].nunique(dropna=False),
            "n_categories_excluding_unknown": sub.loc[
                ~sub["category"].astype(str).str.casefold().eq("unknown"),
                "category",
            ].nunique(dropna=False),
            "n_unique_occurrences": int(sub["n_unique_occurrences"].sum()),
            "n_unknown_occurrences": int(
                sub.loc[
                    sub["category"].astype(str).str.casefold().eq("unknown"),
                    "n_unique_occurrences",
                ].sum()
            ),
        })
    return pd.DataFrame(summary_rows), breakdown


def collapse_for_plot(sub: pd.DataFrame, top_n: int) -> pd.DataFrame:
    top_n = max(1, int(top_n))
    known = sub.loc[~sub["category"].astype(str).str.casefold().eq("unknown")].copy()
    unknown = sub.loc[sub["category"].astype(str).str.casefold().eq("unknown")].copy()
    top = known.head(top_n).copy()
    others_count = int(known.iloc[top_n:]["n_unique_occurrences"].sum())
    pieces = [top]
    if others_count:
        pieces.append(pd.DataFrame([{
            "dimension": sub["dimension"].iloc[0],
            "category": "Others",
            "n_unique_occurrences": others_count,
        }]))
    if not unknown.empty:
        pieces.append(pd.DataFrame([{
            "dimension": sub["dimension"].iloc[0],
            "category": "Unknown",
            "n_unique_occurrences": int(unknown["n_unique_occurrences"].sum()),
        }]))
    out = pd.concat(pieces, ignore_index=True)
    out["percent_unique_occurrences"] = (
        out["n_unique_occurrences"] / out["n_unique_occurrences"].sum() * 100
    )
    return out


def plot_donut(sub: pd.DataFrame, output_path: Path, top_n: int, font_scale: float) -> None:
    plot_data = collapse_for_plot(sub, top_n)
    labels = plot_data["category"].astype(str).tolist()
    values = plot_data["n_unique_occurrences"].astype(int).tolist()
    colors = plt.cm.tab20(np.linspace(0, 1, max(3, len(values))))
    for i, label in enumerate(labels):
        if label in {"Others", "Unknown"}:
            colors[i] = (0.62, 0.62, 0.62, 1.0) if label == "Others" else (0.82, 0.82, 0.82, 1.0)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    wedges, _ = ax.pie(
        values,
        startangle=90,
        colors=colors,
        wedgeprops={"width": 0.42, "edgecolor": "white", "linewidth": 0.8},
    )
    ax.set(aspect="equal")
    total = sum(values)
    dimension = sub["dimension"].iloc[0]
    ax.set_title(f"{dimension}: unique occurrence breakdown", fontsize=11 * font_scale)
    ax.text(
        0,
        0,
        f"n = {total:,}",
        ha="center",
        va="center",
        fontsize=10 * font_scale,
    )
    legend_labels = [
        f"{label} ({count:,}; {count / total * 100:.1f}%)"
        for label, count in zip(labels, values)
    ]
    ax.legend(
        wedges,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=7.5 * font_scale,
    )
    fig.tight_layout()
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_top_sh_bar(df: pd.DataFrame, output_path: Path, table_path: Path, top_n: int, font_scale: float) -> None:
    """Write a ranked SH occurrence table and a horizontal bar chart."""
    top_n = max(1, int(top_n))
    counts = (
        df["sh_id"]
        .map(clean_label)
        .value_counts(dropna=False)
        .rename_axis("sh_id")
        .reset_index(name="n_unique_occurrences")
    )
    counts["rank"] = np.arange(1, len(counts) + 1)
    counts["percent_unique_occurrences"] = counts["n_unique_occurrences"] / len(df) * 100
    counts.to_csv(table_path, sep="\t", index=False)

    plot_data = counts.head(top_n).iloc[::-1].copy()
    height = max(6.0, 0.18 * len(plot_data) + 1.8)
    fig, ax = plt.subplots(figsize=(8.0, height))
    color = "#8b2f67"
    ax.barh(plot_data["sh_id"], plot_data["n_unique_occurrences"], color=color, alpha=0.86)
    ax.set_xlabel("Unique occurrences", fontsize=10 * font_scale)
    ax.set_ylabel("Fungal SH", fontsize=10 * font_scale)
    ax.set_title(
        f"Top {min(top_n, len(counts))} fungal SHs by unique occurrence count",
        fontsize=11 * font_scale,
    )
    ax.tick_params(axis="both", labelsize=7.5 * font_scale)
    ax.grid(axis="x", color="#bdbdbd", linewidth=0.45, alpha=0.7)
    ax.set_axisbelow(True)
    xmax = plot_data["n_unique_occurrences"].max()
    for y, value in enumerate(plot_data["n_unique_occurrences"]):
        ax.text(
            value + xmax * 0.01,
            y,
            f"{int(value):,}",
            va="center",
            ha="left",
            fontsize=7.0 * font_scale,
        )
    ax.set_xlim(0, xmax * 1.14)
    fig.tight_layout()
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def write_run_files(args: argparse.Namespace, paths: dict[str, Path], start: datetime, end: datetime) -> None:
    output_dir = paths["output_dir"]
    command = " ".join(shlex.quote(part) for part in sys.argv)
    (output_dir / "run_command.sh").write_text("#!/usr/bin/env bash\n" + command + "\n")
    pd.DataFrame([
        {"key": "started_at_utc", "value": start.isoformat()},
        {"key": "finished_at_utc", "value": end.isoformat()},
        {"key": "elapsed_seconds", "value": f"{(end - start).total_seconds():.3f}"},
        {"key": "python", "value": sys.version.replace("\n", " ")},
        {"key": "python_executable", "value": sys.executable},
        {"key": "platform", "value": platform.platform()},
        {"key": "input", "value": rel(paths["input"])},
        {"key": "specificity_result_dir", "value": rel(paths["result_dir"])},
        {"key": "plant_taxonomy", "value": rel(paths["taxonomy"])},
        {"key": "output_dir", "value": rel(output_dir)},
        {"key": "plot_top_n", "value": args.plot_top_n},
        {"key": "top_sh_bar_n", "value": args.top_sh_bar_n},
    ]).to_csv(output_dir / "run_summary.tsv", sep="\t", index=False)


def main() -> int:
    start = datetime.now(timezone.utc)
    args = parse_args()
    if not math.isfinite(args.figure_font_scale) or args.figure_font_scale <= 0:
        raise ValueError("--figure-font-scale must be a positive finite value.")
    paths = resolve_paths(args)
    require_files(paths)
    output_dir = paths["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    cleaned, removal_summary = clean_occurrences(paths)
    summary, breakdown = build_breakdowns(cleaned)

    stem = paths["stem"].name
    cleaned_path = output_dir / f"{stem}_cleaned_for_data_property.tsv"
    cleaned.to_csv(cleaned_path, sep="\t", index=False)
    removal_summary.to_csv(output_dir / "cleaning_summary.tsv", sep="\t", index=False)
    summary.to_csv(output_dir / "category_count_summary.tsv", sep="\t", index=False)
    breakdown.to_csv(output_dir / "unique_occurrence_breakdown_by_category.tsv", sep="\t", index=False)

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    plot_top_sh_bar(
        cleaned,
        figure_dir / f"top_{args.top_sh_bar_n}_fungal_sh_unique_occurrence_barplot.pdf",
        output_dir / "fungal_sh_unique_occurrence_ranking.tsv",
        args.top_sh_bar_n,
        args.figure_font_scale,
    )
    with PdfPages(figure_dir / "unique_occurrence_breakdown_all_dimensions.pdf") as pdf:
        for dimension, sub in breakdown.groupby("dimension", sort=False):
            path = figure_dir / f"unique_occurrence_breakdown_{dimension}.pdf"
            plot_donut(sub, path, args.plot_top_n, args.figure_font_scale)
            plot_data = collapse_for_plot(sub, args.plot_top_n)
            fig, ax = plt.subplots(figsize=(6.2, 5.2))
            wedges, _ = ax.pie(
                plot_data["n_unique_occurrences"],
                startangle=90,
                wedgeprops={"width": 0.42, "edgecolor": "white", "linewidth": 0.8},
            )
            ax.set(aspect="equal")
            ax.set_title(f"{dimension}: unique occurrence breakdown", fontsize=11 * args.figure_font_scale)
            ax.text(
                0,
                0,
                f"n = {int(plot_data['n_unique_occurrences'].sum()):,}",
                ha="center",
                va="center",
                fontsize=10 * args.figure_font_scale,
            )
            ax.legend(
                wedges,
                [
                    f"{row.category} ({int(row.n_unique_occurrences):,}; {row.percent_unique_occurrences:.1f}%)"
                    for row in plot_data.itertuples(index=False)
                ],
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False,
                fontsize=7.5 * args.figure_font_scale,
            )
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    end = datetime.now(timezone.utc)
    write_run_files(args, paths, start, end)
    print(f"Wrote cleaned table: {cleaned_path}")
    print(f"Cleaned rows: {len(cleaned):,}; SH: {cleaned['sh_id'].nunique():,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
