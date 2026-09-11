#!/usr/bin/env python3
"""Create overview maps and continent-composition plots for GlobalFungi/Helotiales data.

The script summarizes four nested data definitions:

1. All GlobalFungi occurrence records.
2. All Helotiales occurrence records.
3. All Helotiales occurrence records from root samples.
4. Reliable root-host unique occurrences after data-property cleaning.

Large GlobalFungi-derived tables are read in chunks so that the overview can be
regenerated without loading the full occurrence table into memory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shlex
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

CACHE_DIR = Path(tempfile.gettempdir()) / "helotiales_data_overview_matplotlib"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_DIR))

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
UNKNOWN_TOKENS = {"", "NA", "N/A", "NAN", "NONE", "NULL", "UNKNOWN", "<NA>"}
CONTINENT_ORDER = [
    "Asia",
    "Europe",
    "North America",
    "South America",
    "Africa",
    "Australia",
    "Antarctica",
    "Pacific Ocean",
    "Atlantic Ocean",
    "Arctic Ocean",
    "Indian Ocean",
    "Unknown",
    "Other",
]
CONTINENT_COLORS = {
    "Asia": "#7F3C8D",
    "Europe": "#11A579",
    "North America": "#3969AC",
    "South America": "#F2B701",
    "Africa": "#E73F74",
    "Australia": "#80BA5A",
    "Antarctica": "#E68310",
    "Pacific Ocean": "#008695",
    "Atlantic Ocean": "#CF1C90",
    "Arctic Ocean": "#F97B72",
    "Indian Ocean": "#A5AA99",
    "Unknown": "#999999",
    "Other": "#CCCCCC",
}


DATASETS = [
    {
        "key": "all_globalfungi_occurrences",
        "label": "All GlobalFungi occurrences",
        "unit_label": "occurrences",
    },
    {
        "key": "all_helotiales_occurrences",
        "label": "All Helotiales occurrences",
        "unit_label": "occurrences",
    },
    {
        "key": "root_helotiales_occurrences",
        "label": "Root-sample Helotiales occurrences",
        "unit_label": "occurrences",
    },
    {
        "key": "reliable_root_host_unique_occurrences",
        "label": "Reliable root-host unique occurrences",
        "unit_label": "unique occurrences",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create overview maps and continent proportion plots for Helotiales data."
    )
    parser.add_argument(
        "--globalfungi-results",
        type=Path,
        default=PROJECT_DIR / "1_GlobalFungi" / "results_helotiales",
        help="Directory containing 1_GlobalFungi outputs.",
    )
    parser.add_argument(
        "--globalfungi-data",
        type=Path,
        default=PROJECT_DIR / "1_GlobalFungi" / "data",
        help="Directory containing raw GlobalFungi input files.",
    )
    parser.add_argument(
        "--globalfungi-abundance",
        type=Path,
        default=None,
        help=(
            "GlobalFungi SH abundance matrix. Default: "
            "<globalfungi-data>/GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz."
        ),
    )
    parser.add_argument(
        "--globalfungi-metadata",
        type=Path,
        default=None,
        help=(
            "GlobalFungi sample metadata table. Default: "
            "<globalfungi-data>/GlobalFungi_5_sample_metadata.txt.gz."
        ),
    )
    parser.add_argument(
        "--all-occurrences",
        type=Path,
        default=None,
        help="All Helotiales occurrence table. Default: <globalfungi-results>/helotiales_occurrences_merged.tsv.",
    )
    parser.add_argument(
        "--reliable-host-occurrences",
        type=Path,
        default=None,
        help=(
            "Reliable root-host unique-occurrence table. By default the script first "
            "uses the 0_Data_Property cleaned table if present, then falls back to "
            "<globalfungi-results>/helotiales_root_occurrences_with_hosts.tsv."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR / "outputs",
        help="Output directory.",
    )
    parser.add_argument(
        "--world-geojson",
        type=Path,
        default=SCRIPT_DIR / "assets" / "ne_110m_admin_0_countries.geojson",
        help="Natural Earth country-boundary GeoJSON.",
    )
    parser.add_argument("--chunksize", type=int, default=200_000)
    parser.add_argument("--figure-font-scale", type=float, default=1.5)
    parser.add_argument("--point-size", type=float, default=10.0)
    parser.add_argument("--point-alpha", type=float, default=0.55)
    args = parser.parse_args()
    if args.chunksize < 1:
        parser.error("--chunksize must be positive")
    if args.figure_font_scale <= 0:
        parser.error("--figure-font-scale must be positive")
    if args.point_size <= 0:
        parser.error("--point-size must be positive")
    if not 0 < args.point_alpha <= 1:
        parser.error("--point-alpha must be in (0, 1]")
    return args


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_DIR.resolve()).as_posix()
    except Exception:
        return path.name


def normalize_text(value: object) -> str:
    if pd.isna(value):
        return "Unknown"
    text = str(value).strip()
    if text.upper() in UNKNOWN_TOKENS:
        return "Unknown"
    return text


def valid_coordinate_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["latitude"] = pd.to_numeric(out["latitude"], errors="coerce")
    out["longitude"] = pd.to_numeric(out["longitude"], errors="coerce")
    return out.loc[
        out["latitude"].between(-90, 90) & out["longitude"].between(-180, 180),
        ["latitude", "longitude"],
    ]


def update_counts(
    store: dict[str, object],
    chunk: pd.DataFrame,
    unit_increment: int | None = None,
) -> None:
    n_rows = len(chunk) if unit_increment is None else int(unit_increment)
    store["n_records"] = int(store["n_records"]) + n_rows
    coords = valid_coordinate_frame(chunk)
    store["n_records_with_coordinates"] = int(store["n_records_with_coordinates"]) + len(coords)
    if not coords.empty:
        store["coordinate_pairs"].update(map(tuple, coords.drop_duplicates().to_numpy()))
    continents = chunk["continent"].map(normalize_text).fillna("Unknown")
    store["continent_counts"].update(continents)


def new_store() -> dict[str, object]:
    return {
        "n_records": 0,
        "n_records_with_coordinates": 0,
        "coordinate_pairs": set(),
        "continent_counts": Counter(),
    }


def read_globalfungi_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(
        path,
        sep="\t",
        usecols=["sample_ID", "latitude", "longitude", "continent"],
        dtype=str,
        encoding="utf-8-sig",
    )
    metadata = metadata.drop_duplicates("sample_ID", keep="first").set_index("sample_ID")
    metadata["latitude"] = pd.to_numeric(metadata["latitude"], errors="coerce")
    metadata["longitude"] = pd.to_numeric(metadata["longitude"], errors="coerce")
    metadata["continent"] = metadata["continent"].map(normalize_text)
    return metadata


def read_globalfungi_all_occurrences(abundance_path: Path, metadata_path: Path) -> dict[str, object]:
    store = new_store()
    metadata = read_globalfungi_metadata(metadata_path)
    sample_counts: dict[str, int] = {}
    opener = gzip.open if abundance_path.suffix == ".gz" else open
    with opener(abundance_path, "rb") as handle:
        header = handle.readline()
        if not header.startswith(b"sample_ID\t"):
            raise ValueError(f"Unexpected GlobalFungi abundance header in {abundance_path}")
        n_sh_columns = header.count(b"\t")
        for raw in handle:
            tab = raw.find(b"\t")
            if tab < 1:
                continue
            sample_id = raw[:tab].decode("utf-8", "replace")
            row_data = raw[tab:].rstrip(b"\r\n") + b"\t"
            n_positive = n_sh_columns - row_data.count(b"\t0")
            if n_positive > 0:
                sample_counts[sample_id] = n_positive
    if not sample_counts:
        return store

    counts = pd.Series(sample_counts, name="occurrence_count", dtype="int64")
    merged = metadata.join(counts, how="inner")
    store["n_records"] = int(merged["occurrence_count"].sum())
    valid = merged.loc[
        merged["latitude"].between(-90, 90) & merged["longitude"].between(-180, 180)
    ].copy()
    store["n_records_with_coordinates"] = int(valid["occurrence_count"].sum())
    if not valid.empty:
        store["coordinate_pairs"].update(
            map(tuple, valid[["latitude", "longitude"]].drop_duplicates().to_numpy())
        )
    for continent, total in merged.groupby("continent")["occurrence_count"].sum().items():
        store["continent_counts"][continent] += int(total)
    store["n_samples_with_occurrences"] = int(len(merged))
    store["n_samples_with_coordinates"] = int(len(valid))
    return store


def read_all_and_root_occurrences(path: Path, chunksize: int) -> tuple[dict[str, object], dict[str, object]]:
    all_store = new_store()
    root_store = new_store()
    required = ["sample_type", "latitude", "longitude", "continent"]
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=required,
        dtype=str,
        chunksize=chunksize,
        encoding="utf-8-sig",
        quoting=csv.QUOTE_NONE,
        on_bad_lines="skip",
        low_memory=False,
    ):
        update_counts(all_store, chunk)
        sample_type = chunk["sample_type"].astype("string").str.strip().str.casefold()
        root = chunk.loc[sample_type.eq("root")].copy()
        update_counts(root_store, root)
    return all_store, root_store


def read_reliable_host_occurrences(path: Path) -> dict[str, object]:
    store = new_store()
    cols = pd.read_csv(path, sep="\t", nrows=0, encoding="utf-8-sig").columns
    wanted = [column for column in ["latitude", "longitude", "continent", "unique_occurrence_count"] if column in cols]
    df = pd.read_csv(path, sep="\t", usecols=wanted, dtype=str, encoding="utf-8-sig")
    if "unique_occurrence_count" in df.columns:
        weights = pd.to_numeric(df["unique_occurrence_count"], errors="coerce").fillna(1).astype(int)
        expanded_continents = Counter()
        for continent, weight in zip(df["continent"].map(normalize_text), weights):
            expanded_continents[continent] += int(weight)
        store["n_records"] = int(weights.sum())
        coords = valid_coordinate_frame(df)
        store["n_records_with_coordinates"] = int(weights.loc[coords.index].sum()) if not coords.empty else 0
        store["coordinate_pairs"].update(map(tuple, coords.drop_duplicates().to_numpy()))
        store["continent_counts"].update(expanded_continents)
    else:
        update_counts(store, df)
    return store


def draw_world_boundaries(ax: plt.Axes, geojson_path: Path) -> int:
    with geojson_path.open("r", encoding="utf-8") as handle:
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
                breaks = np.flatnonzero(np.abs(np.diff(points[:, 0])) > 180) + 1
                for part in np.split(points[:, :2], breaks):
                    if len(part) >= 2:
                        segments.append(part)
    ax.add_collection(LineCollection(segments, colors="#686868", linewidths=0.45, alpha=0.9, zorder=1))
    return len(segments)


def apply_font_scale(fig: plt.Figure, scale: float) -> None:
    for text in fig.findobj(match=plt.Text):
        text.set_fontsize(text.get_fontsize() * scale)


def plot_map(store: dict[str, object], dataset: dict[str, str], output: Path, geojson_path: Path, font_scale: float, point_size: float, point_alpha: float) -> None:
    coords = pd.DataFrame(list(store["coordinate_pairs"]), columns=["latitude", "longitude"])
    fig, ax = plt.subplots(figsize=(11.5, 6.0))
    ax.set_facecolor("#f7fbfc")
    draw_world_boundaries(ax, geojson_path)
    if not coords.empty:
        ax.scatter(
            coords["longitude"],
            coords["latitude"],
            s=point_size,
            alpha=point_alpha,
            color="#8B2F67",
            edgecolors="none",
            zorder=2,
        )
    ax.set(xlim=(-180, 180), ylim=(-90, 90), xlabel="Longitude", ylabel="Latitude")
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.grid(color="#aaaaaa", linewidth=0.35, alpha=0.45)
    ax.set_title(dataset["label"])
    ax.text(
        0.015,
        0.04,
        f"Sampling locations = {len(coords):,}\n{dataset['unit_label']} = {int(store['n_records']):,}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.2,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#808080", "linewidth": 0.6, "alpha": 0.88},
        zorder=3,
    )
    apply_font_scale(fig, font_scale)
    fig.tight_layout()
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)


def continent_table(store: dict[str, object], key: str, label: str, unit_label: str) -> pd.DataFrame:
    counts = pd.Series(dict(store["continent_counts"]), dtype="int64")
    unknown_continents = [continent for continent in counts.index if continent not in CONTINENT_ORDER]
    if unknown_continents:
        counts.loc["Other"] = counts.get("Other", 0) + counts.loc[unknown_continents].sum()
        counts = counts.drop(index=unknown_continents)
    ordered = [continent for continent in CONTINENT_ORDER if continent in counts.index]
    counts = counts.loc[ordered]
    total = int(counts.sum())
    return pd.DataFrame({
        "dataset": key,
        "dataset_label": label,
        "unit": unit_label,
        "continent": counts.index,
        "count": counts.values,
        "proportion": counts.values / total if total else np.nan,
    })


def plot_continent_donut(table: pd.DataFrame, dataset: dict[str, str], output: Path, font_scale: float) -> None:
    data = table.loc[table["dataset"].eq(dataset["key"])].copy()
    order = {continent: i for i, continent in enumerate(CONTINENT_ORDER)}
    data["continent_order"] = data["continent"].map(order).fillna(len(CONTINENT_ORDER)).astype(int)
    data = data.sort_values("continent_order")
    colors = [CONTINENT_COLORS.get(continent, CONTINENT_COLORS["Other"]) for continent in data["continent"]]
    labels = [f"{row.continent}\n{row.proportion:.1%}" for row in data.itertuples()]
    fig, ax = plt.subplots(figsize=(7.1, 6.4))
    wedges, _ = ax.pie(
        data["count"],
        labels=None,
        startangle=90,
        counterclock=False,
        colors=colors,
        wedgeprops={"width": 0.38, "edgecolor": "white", "linewidth": 0.8},
    )
    outside_labels: list[dict[str, object]] = []
    for wedge, label, proportion in zip(wedges, labels, data["proportion"]):
        angle = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)
        if proportion >= 0.045:
            radius = 0.78
            ax.text(
                radius * np.cos(angle),
                radius * np.sin(angle),
                label,
                ha="center",
                va="center",
                fontsize=7.0,
                color="#111111",
            )
        else:
            side = 1 if np.cos(angle) >= 0 else -1
            outside_labels.append({
                "label": label,
                "angle": angle,
                "side": side,
                "target": (0.93 * np.cos(angle), 0.93 * np.sin(angle)),
            })

    placed_outside = []
    for side in (-1, 1):
        items = sorted(
            [item for item in outside_labels if item["side"] == side],
            key=lambda item: float(item["angle"]),
        )
        previous_angle = None
        radial_level = 0
        for item in items:
            angle = float(item["angle"])
            if previous_angle is None or angle - previous_angle > 0.18:
                radial_level = 0
            else:
                radial_level += 1
            item["label_radius"] = 1.10 + 0.105 * radial_level
            placed_outside.append(item)
            previous_angle = angle

    for item in placed_outside:
        side = int(item["side"])
        angle = float(item["angle"])
        label_radius = float(item["label_radius"])
        line_end_radius = max(1.02, label_radius - 0.055)
        x = label_radius * np.cos(angle)
        y = label_radius * np.sin(angle)
        target_x, target_y = item["target"]
        ax.plot(
            [target_x, line_end_radius * np.cos(angle)],
            [target_y, line_end_radius * np.sin(angle)],
            color="#666666",
            lw=0.36,
            solid_capstyle="round",
            zorder=3,
        )
        ax.text(
            x,
            y,
            str(item["label"]),
            ha="left" if side > 0 else "right",
            va="center",
            fontsize=4.5,
            linespacing=0.92,
            zorder=4,
        )
    ax.text(0, 0, f"n = {int(data['count'].sum()):,}", ha="center", va="center", fontsize=10)
    ax.set_title(dataset["label"], pad=20)
    ax.set(aspect="equal")
    ax.set_xlim(-1.86, 1.86)
    ax.set_ylim(-1.48, 2.18)
    apply_font_scale(fig, font_scale)
    fig.tight_layout()
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)


def write_runtime_versions(output: Path) -> None:
    rows = [
        {"category": "python", "name": "Python", "version": platform.python_version()},
        {"category": "python_package", "name": "pandas", "version": pd.__version__},
        {"category": "python_package", "name": "numpy", "version": np.__version__},
        {"category": "python_package", "name": "matplotlib", "version": matplotlib.__version__},
    ]
    for dist in importlib_metadata.distributions():
        name = dist.metadata.get("Name", "Unknown")
        if name in {"pandas", "numpy", "matplotlib"}:
            continue
        rows.append({"category": "python_package", "name": name, "version": dist.version})
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


def main() -> int:
    started = time.perf_counter()
    started_at = datetime.now(timezone.utc).astimezone()
    args = parse_args()
    args.globalfungi_results = args.globalfungi_results.resolve()
    args.globalfungi_data = args.globalfungi_data.resolve()
    globalfungi_abundance = (
        args.globalfungi_abundance
        or args.globalfungi_data / "GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz"
    ).resolve()
    globalfungi_metadata = (
        args.globalfungi_metadata
        or args.globalfungi_data / "GlobalFungi_5_sample_metadata.txt.gz"
    ).resolve()
    all_path = (args.all_occurrences or args.globalfungi_results / "helotiales_occurrences_merged.tsv").resolve()
    default_cleaned = PROJECT_DIR / "0_Data_Property" / "outputs" / "sh150_family30_10000_lat20_minblock100" / "helotiales_root_occurrences_with_hosts_cleaned_for_data_property.tsv"
    reliable_path = (args.reliable_host_occurrences or (default_cleaned if default_cleaned.is_file() else args.globalfungi_results / "helotiales_root_occurrences_with_hosts.tsv")).resolve()
    args.output_dir = args.output_dir.resolve()
    args.world_geojson = args.world_geojson.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = args.output_dir / "figures"
    figures.mkdir(exist_ok=True)
    log = args.output_dir / "workflow.log"
    command = shlex.join([sys.executable, *sys.argv])
    log.write_text(f"Started: {started_at.isoformat()}\nCommand: {command}\n", encoding="utf-8")

    for path, label in [
        (globalfungi_abundance, "GlobalFungi abundance matrix"),
        (globalfungi_metadata, "GlobalFungi sample metadata"),
        (all_path, "all Helotiales occurrences"),
        (reliable_path, "reliable host occurrences"),
        (args.world_geojson, "world GeoJSON"),
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")

    globalfungi_store = read_globalfungi_all_occurrences(globalfungi_abundance, globalfungi_metadata)
    all_store, root_store = read_all_and_root_occurrences(all_path, args.chunksize)
    reliable_store = read_reliable_host_occurrences(reliable_path)
    stores = {
        "all_globalfungi_occurrences": globalfungi_store,
        "all_helotiales_occurrences": all_store,
        "root_helotiales_occurrences": root_store,
        "reliable_root_host_unique_occurrences": reliable_store,
    }

    summary_rows = []
    continent_tables = []
    for dataset in DATASETS:
        store = stores[dataset["key"]]
        summary_rows.append({
            "dataset": dataset["key"],
            "dataset_label": dataset["label"],
            "unit": dataset["unit_label"],
            "n_units": int(store["n_records"]),
            "n_units_with_valid_coordinates": int(store["n_records_with_coordinates"]),
            "n_sampling_locations": len(store["coordinate_pairs"]),
            "n_continents": len([k for k, v in store["continent_counts"].items() if v > 0]),
            "n_samples_with_occurrences": store.get("n_samples_with_occurrences", ""),
            "n_samples_with_coordinates": store.get("n_samples_with_coordinates", ""),
        })
        ctab = continent_table(store, dataset["key"], dataset["label"], dataset["unit_label"])
        continent_tables.append(ctab)
        plot_map(
            store,
            dataset,
            figures / f"{dataset['key']}_world_map.pdf",
            args.world_geojson,
            args.figure_font_scale,
            args.point_size,
            args.point_alpha,
        )
        plot_continent_donut(
            ctab,
            dataset,
            figures / f"{dataset['key']}_continent_proportions.pdf",
            args.figure_font_scale,
        )

    pd.DataFrame(summary_rows).to_csv(args.output_dir / "data_overview_summary.tsv", sep="\t", index=False)
    pd.concat(continent_tables, ignore_index=True).to_csv(args.output_dir / "continent_breakdown.tsv", sep="\t", index=False)
    pd.DataFrame([
        {
            "continent": continent,
            "display_order": i,
            "color": CONTINENT_COLORS[continent],
        }
        for i, continent in enumerate(CONTINENT_ORDER, start=1)
    ]).to_csv(args.output_dir / "continent_color_key.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"setting": "command", "value": command},
        {"setting": "globalfungi_results", "value": rel(args.globalfungi_results)},
        {"setting": "globalfungi_data", "value": rel(args.globalfungi_data)},
        {"setting": "globalfungi_abundance", "value": rel(globalfungi_abundance)},
        {"setting": "globalfungi_metadata", "value": rel(globalfungi_metadata)},
        {"setting": "all_occurrences", "value": rel(all_path)},
        {"setting": "reliable_host_occurrences", "value": rel(reliable_path)},
        {"setting": "world_geojson", "value": rel(args.world_geojson)},
        {"setting": "output_dir", "value": rel(args.output_dir)},
        {"setting": "chunksize", "value": args.chunksize},
        {"setting": "figure_font_scale", "value": args.figure_font_scale},
        {"setting": "point_size", "value": args.point_size},
        {"setting": "point_alpha", "value": args.point_alpha},
    ]).to_csv(args.output_dir / "run_configuration.tsv", sep="\t", index=False)
    (args.output_dir / "run_command.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + command + "\n", encoding="utf-8")
    write_runtime_versions(args.output_dir / "runtime_versions.tsv")
    pd.DataFrame([{
        "start_time": started_at.isoformat(),
        "end_time": datetime.now(timezone.utc).astimezone().isoformat(),
        "elapsed_seconds": time.perf_counter() - started,
        "status": "completed",
        "command": command,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }]).to_csv(args.output_dir / "run_timing.tsv", sep="\t", index=False)
    with log.open("a", encoding="utf-8") as handle:
        handle.write(f"Completed in {time.perf_counter() - started:.3f} seconds\n")
    print(f"Completed data overview: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
