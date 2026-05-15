import json
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA_PATH = Path("../docs/json/temporal_dynamics_disease.json")
OUTPUT_DIR = Path("images/grouped_bar")

DISEASE_TYPES = ["MDRB", "Melanoma", "rCDI"]
STATUS_ORDER = ["persisted", "emerged", "disappeared", "transferred"]
STATUS_LABELS = {
    "persisted": "Persistence",
    "emerged": "Emergence",
    "disappeared": "Disappearance",
    "transferred": "Transfer",
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


def build_dataframe(records: List[Dict]) -> pd.DataFrame:
    df = pd.DataFrame(records).copy()
    if df.empty:
        return df

    df["patients"] = pd.to_numeric(df["patients"], errors="coerce").fillna(0)
    df["status_normalized"] = df["status"].astype(str).str.strip().str.lower()
    df["mgeGroup_normalized"] = df["mgeGroup"].apply(normalize_mge)
    return df


def build_grouped_counts(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df[
            df["status_normalized"].isin(STATUS_ORDER)
            & df["mgeGroup_normalized"].isin(MGE_ORDER)
        ]
        .groupby(["status_normalized", "mgeGroup_normalized"])["patients"]
        .sum()
        .unstack(fill_value=0)
        .reindex(index=STATUS_ORDER, columns=MGE_ORDER, fill_value=0)
    )
    return grouped


def plot_grouped_bar_chart(disease: str, records: List[Dict]) -> None:
    df = build_dataframe(records)
    if df.empty:
        print(f"No data found for {disease}; skipping.")
        return

    grouped_counts = build_grouped_counts(df)

    fig, ax = plt.subplots(figsize=(11, 6.5), dpi=300)
    x = np.arange(len(STATUS_ORDER))
    width = 0.15

    for idx, mge_group in enumerate(MGE_ORDER):
        offset = (idx - (len(MGE_ORDER) - 1) / 2) * width
        y_values = grouped_counts[mge_group].values.astype(float)
        bars = ax.bar(
            x + offset,
            y_values,
            width=width,
            color=MGE_COLORS[mge_group],
            label=mge_group,
            edgecolor="white",
            linewidth=1.0,
        )

        for bar, value in zip(bars, y_values):
            if value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.2,
                f"{int(value)}",
                ha="center",
                va="bottom",
                fontsize=18,
                fontweight="bold",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [STATUS_LABELS[status] for status in STATUS_ORDER],
        fontsize=20,
        fontweight="bold",
    )
    ax.set_xlabel("Transition Type", fontsize=22, fontweight="bold")
    ax.set_ylabel("No. of colocalizations", fontsize=22, fontweight="bold")
    ax.tick_params(axis="y", labelsize=18)
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")

    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(
        title="MGE Class",
        fontsize=20,
        title_fontsize=20,
        loc="upper right",
        frameon=False,
    )

    fig.tight_layout()

    out_prefix = OUTPUT_DIR / f"{disease.lower()}_grouped_bar_transition_by_mge"
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out_prefix.with_suffix('.png')}")
    print(f"Saved {out_prefix.with_suffix('.svg')}")


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Could not find JSON file: {DATA_PATH}")

    data = load_data(DATA_PATH)
    for disease in DISEASE_TYPES:
        plot_grouped_bar_chart(disease, data.get(disease, []))


if __name__ == "__main__":
    main()