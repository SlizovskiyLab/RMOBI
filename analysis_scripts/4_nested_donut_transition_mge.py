import json
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd

DATA_PATH = Path("../docs/json/temporal_dynamics_disease.json")
OUTPUT_DIR = Path("images/nested_donut")

STATUS_ORDER = ["persisted", "emerged", "disappeared", "transferred"]
STATUS_LABELS = {
    "persisted": "Persistence",
    "emerged": "Emergence",
    "disappeared": "Disappearance",
    "transferred": "Transfer",
}
STATUS_COLORS = {
    "persisted": "#BCC8DE",
    "emerged": "#52565a",
    "disappeared": "#7893af",
    "transferred": "#2f3439",
}
MGE_ORDER = ["Plasmids", "Likely IS/TE", "ICE", "Prophages", "Virus"]
MGE_PATTERNS = {
    "Plasmids": ["plasmid"],
    "Likely IS/TE": ["likely is/te", "is/te", "transpos", "insertion sequence"],
    "ICE": ["ice", "integrative conjugative"],
    "Prophages": ["prophage", "proph"],
    "Virus": ["virus", "viral", "vir"],
}
MGE_COLORS = {
    "Plasmids": "#1b5a9d",
    "Likely IS/TE": "#3c7eb7",
    "ICE": "#92beda",
    "Prophages": "#9a2d29",
    "Virus": "#e25e49",
}


def load_data(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def normalize_mge(value: str) -> str:
    value_lower = str(value).strip().lower()
    for group, patterns in MGE_PATTERNS.items():
        if any(pattern in value_lower for pattern in patterns):
            return group
    return "Unclassified"


def build_dataframe(data: Dict) -> pd.DataFrame:
    records: List[Dict] = []
    for disease, disease_records in data.items():
        for record in disease_records:
            record_copy = record.copy()
            record_copy["disease"] = disease
            records.append(record_copy)

    df = pd.DataFrame(records)
    if df.empty:
        return df

    df["patients"] = pd.to_numeric(df["patients"], errors="coerce").fillna(0)
    df["status_normalized"] = df["status"].astype(str).str.strip().str.lower()
    df["mgeGroup_normalized"] = df["mgeGroup"].apply(normalize_mge)
    df = df[
        df["status_normalized"].isin(STATUS_ORDER)
        & df["mgeGroup_normalized"].isin(MGE_ORDER)
    ].copy()
    return df


def autopct_threshold(threshold: float = 3.0):
    def _autopct(pct: float) -> str:
        return f"{pct:.1f}%" if pct >= threshold else ""

    return _autopct


def build_nested_rings(df: pd.DataFrame):
    inner_series = (
        df.groupby("status_normalized")["patients"]
        .sum()
        .reindex(STATUS_ORDER, fill_value=0)
    )
    inner_series = inner_series[inner_series > 0]

    outer_values = []
    outer_labels = []
    outer_colors = []

    for status in inner_series.index:
        status_df = df[df["status_normalized"] == status]
        mge_series = (
            status_df.groupby("mgeGroup_normalized")["patients"]
            .sum()
            .reindex(MGE_ORDER, fill_value=0)
        )
        for mge_group, value in mge_series.items():
            if value <= 0:
                continue
            outer_values.append(value)
            outer_labels.append(mge_group)
            outer_colors.append(MGE_COLORS[mge_group])

    return inner_series, outer_values, outer_labels, outer_colors


def plot_nested_donut(df: pd.DataFrame) -> None:
    if df.empty:
        print("No valid records found for nested donut plot; skipping.")
        return

    inner_series, outer_values, outer_labels, outer_colors = build_nested_rings(df)
    total_patients = int(df["patients"].sum())
    total_colocalizations = len(df)

    fig, ax = plt.subplots(figsize=(14, 10), dpi=300)

    inner_wedges, _, inner_autotexts = ax.pie(
        inner_series.values,
        radius=1.0,
        labels=None,
        startangle=90,
        autopct=autopct_threshold(4.0),
        pctdistance=0.86,
        colors=[STATUS_COLORS[status] for status in inner_series.index],
        textprops={"fontsize": 20, "weight": "bold"},
        wedgeprops={"width": 0.30, "edgecolor": "white", "linewidth": 1.5},
    )

    outer_wedges, _, outer_autotexts = ax.pie(
        outer_values,
        radius=1.32,
        labels=None,
        startangle=90,
        autopct=autopct_threshold(3.0),
        pctdistance=0.88,
        colors=outer_colors,
        textprops={"fontsize": 20, "weight": "bold"},
        wedgeprops={"width": 0.28, "edgecolor": "white", "linewidth": 1.4},
    )

    for autotext in list(inner_autotexts) + list(outer_autotexts):
        autotext.set_color("white")
        autotext.set_weight("bold")

    transition_lines = "\n".join(
        f"{STATUS_LABELS[s]}: {int(inner_series[s])}"
        for s in inner_series.index
    )
    center_text = f"Total \n{transition_lines}"
    ax.text(
        0,
        0,
        center_text,
        ha="center",
        va="center",
        fontsize=24,
        fontweight="bold",
        linespacing=1.6,
    )

    status_handles = [
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=STATUS_COLORS[status], markersize=14)
        for status in inner_series.index
    ]
    status_labels = [STATUS_LABELS[status] for status in inner_series.index]

    mge_handles = [
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=MGE_COLORS[mge], markersize=14)
        for mge in MGE_ORDER
    ]

    legend_status = ax.legend(
        status_handles,
        status_labels,
        title="Transition Type",
        loc="center left",
        bbox_to_anchor=(0.98, 0.78),
        frameon=False,
        fontsize=18,
        title_fontsize=18,
    )
    ax.add_artist(legend_status)

    ax.legend(
        mge_handles,
        MGE_ORDER,
        title="MGE Class",
        loc="center left",
        bbox_to_anchor=(0.98, 0.28),
        frameon=False,
        fontsize=18,
        title_fontsize=18,
    )

    ax.set(aspect="equal")
    plt.tight_layout()

    out_prefix = OUTPUT_DIR / "nested_donut_transition_mge"
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out_prefix.with_suffix('.png')}")
    print(f"Saved {out_prefix.with_suffix('.svg')}")
    print(f"Total colocalizations plotted: {total_colocalizations}")
    print(f"Total patient-weighted value plotted: {total_patients}")


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Could not find JSON file: {DATA_PATH}")

    data = load_data(DATA_PATH)
    df = build_dataframe(data)
    plot_nested_donut(df)


if __name__ == "__main__":
    main()