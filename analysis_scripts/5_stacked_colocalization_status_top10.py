import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

DATA_PATH = PROJECT_ROOT / "docs/json/temporal_dynamics_disease.json"
CLASSIFICATION_PATH = PROJECT_ROOT / "data/MGE_total_classification.csv"
OUTPUT_DIR = PROJECT_ROOT / "images/stacked"
SUMMARY_OUTPUT_PATH = (
    PROJECT_ROOT
    / "docs/output/top_entities/most_prominent_colocalizations_by_disease.csv"
)

STATUS_ORDER = ["emerged", "disappeared", "persisted", "transferred"]
STATUS_LABELS = {
    "emerged": "Emergence",
    "disappeared": "Disappearance",
    "persisted": "Persistence",
    "transferred": "Transfer",
}

STATUS_COLORS = {
    "emerged": "#9a2d29",    # 
    "disappeared": "#1b5a9d",  # 
    "persisted": "#92beda",    # 
    "transferred": "#e25e49",  # 
}

# Handles common misspellings/variants and maps to canonical status labels.
STATUS_ALIASES = {
    "emerged": "emerged",
    "emergence": "emerged",
    "disappeared": "disappeared",
    "disappearance": "disappeared",
    "dissappeared": "disappeared",
    "persisted": "persisted",
    "persistence": "persisted",
    "persistance": "persisted",
    "transferred": "transferred",
    "transfer": "transferred",
}


def normalize_status(value: str) -> str:
    key = str(value).strip().lower()
    return STATUS_ALIASES.get(key, key)


def shorten(text: str, max_len: int = 42) -> str:
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def _clean_field(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _normalize_classification(value: object) -> str:
    return re.sub(r"\s+", "", _clean_field(value).lower())


def _extract_ice_name(identifier: str) -> str:
    parts = [p.strip() for p in identifier.split("|")]
    if len(parts) >= 3 and parts[2]:
        return parts[2]
    return identifier


def build_mge_name_lookup(csv_path: Path) -> Dict[str, str]:
    if not csv_path.exists():
        raise FileNotFoundError(f"Could not find classification CSV file: {csv_path}")

    df = pd.read_csv(csv_path, encoding="utf-8", low_memory=False)
    required = {"IDs", "final_classification", "socus", "sig_seq"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"Missing required columns in classification CSV: {missing}")

    lookup: Dict[str, str] = {}
    for _, row in df.iterrows():
        mge_id = _clean_field(row["IDs"])
        if not mge_id:
            continue

        classification = _normalize_classification(row["final_classification"])
        socus = _clean_field(row["socus"])
        sig_seq = _clean_field(row["sig_seq"])

        label = ""
        # Rule 2: likely IS/TE should be named by sig_seq for plasmids/prophages.
        if classification == "likelyis/te":
            label = sig_seq or mge_id
        # Rule 1: plasmid/prophage classification named using socus.
        elif classification == "plasmid":
            label = socus or mge_id
        elif classification == "prophage":
            label = socus or mge_id
        # Rule 3: ICE/ICEberg named by the 3rd pipe-delimited token in IDs.
        elif classification in {"ice", "iceberg"}:
            label = _extract_ice_name(mge_id)
        # Rule 4: virus classified as gene:vir.
        elif classification == "virus":
            label = "gene:vir"
        # Rule 5: keep full first-column IDs for Inc/Colicin plasmid and replicon.
        elif classification in {"inc_plasmid", "colicin_plasmid", "replicon"}:
            label = mge_id

        if label:
            lookup[mge_id] = label

    return lookup


def format_colocalization_label(coloc: str, mge_name_lookup: Dict[str, str]) -> str:
    text = str(coloc).strip()
    pieces = re.split(r"[\u2013-]", text, maxsplit=1)
    if len(pieces) != 2:
        return shorten(text, 34)

    arg_name = pieces[0].strip()
    rhs = pieces[1].strip()
    rhs_parts = rhs.split(":")
    mapped = mge_name_lookup.get(rhs, "")
    if mapped:
        return shorten(f"{arg_name}-{mapped}", 34)

    if len(rhs_parts) < 3:
        return shorten(f"{arg_name}-{rhs}", 34)

    mge_type = rhs_parts[1].strip().lower()
    number = rhs_parts[2].strip()

    if mge_type == "plasmid":
        return f"{arg_name}-plas:{number}"
    if mge_type in {"proph", "prophage"}:
        return f"{arg_name}-prop:{number}"
    return shorten(f"{arg_name}-{rhs}", 34)


def build_top10_table(records: List[Dict]) -> Tuple[pd.DataFrame, List[str]]:
    df = pd.DataFrame(records)
    required = {"colocalization", "status", "patients"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")

    df = df.copy()
    df["status_norm"] = df["status"].apply(normalize_status)
    df = df[df["status_norm"].isin(STATUS_ORDER)].copy()
    df["patients"] = pd.to_numeric(df["patients"], errors="coerce").fillna(0)

    grouped = (
        df.groupby(["colocalization", "status_norm"], as_index=False)["patients"]
        .sum()
        .sort_values("patients", ascending=False)
    )

    totals = (
        grouped.groupby("colocalization", as_index=True)["patients"]
        .sum()
        .sort_values(ascending=False)
    )
    top10 = totals.head(10).index.tolist()

    pivot = grouped[grouped["colocalization"].isin(top10)].pivot_table(
        index="status_norm",
        columns="colocalization",
        values="patients",
        aggfunc="sum",
        fill_value=0,
    )

    pivot = pivot.reindex(STATUS_ORDER, fill_value=0)
    pivot = pivot.reindex(columns=top10, fill_value=0)
    return pivot, top10


def build_prominent_colocalizations_table(
    data: Dict[str, List[Dict]], top_n: int = 10
) -> pd.DataFrame:
    """Rank ARG-MGE pairs by recipient patient counts across disease cohorts.

    Only recognized temporal-transition records are included. This explicitly
    excludes donor-only observations while retaining transferred observations,
    which were detected in recipients after FMT.
    """
    diseases = list(data.keys())
    aggregated: Dict[str, Dict[str, object]] = {}

    for disease, records in data.items():
        for record in records:
            status = normalize_status(record.get("status", ""))
            if status not in STATUS_ORDER:
                continue

            colocalization = _clean_field(record.get("colocalization", ""))
            if not colocalization:
                continue

            entry = aggregated.setdefault(
                colocalization,
                {
                    "patients": {name: 0.0 for name in diseases},
                    "studies": {name: set() for name in diseases},
                },
            )
            patient_counts = entry["patients"]
            study_sets = entry["studies"]

            patients = pd.to_numeric(record.get("patients", 0), errors="coerce")
            if pd.isna(patients):
                patients = 0
            patient_counts[disease] += float(patients)

            disease_studies = record.get("diseaseStudyCounts", {}).get(disease, {})
            if isinstance(disease_studies, dict):
                study_sets[disease].update(
                    str(study) for study in disease_studies if _clean_field(study)
                )
            else:
                study_counts = record.get("studyCounts", {})
                if isinstance(study_counts, dict):
                    study_sets[disease].update(
                        str(study) for study in study_counts if _clean_field(study)
                    )
                else:
                    study_sets[disease].update(
                        str(study)
                        for study in record.get("studies", [])
                        if _clean_field(study)
                    )

    rows = []
    for colocalization, entry in aggregated.items():
        patient_counts = entry["patients"]
        study_sets = entry["studies"]
        all_studies = set().union(*(study_sets[disease] for disease in diseases))
        row = {
            "Colocalization": colocalization,
            "Total patients": sum(patient_counts.values()),
            "Total source studies": len(all_studies),
        }
        for disease in diseases:
            row[f"{disease} patients"] = patient_counts[disease]
            row[f"{disease} source studies"] = len(study_sets[disease])
        rows.append(row)

    columns = [
        "Rank",
        "Colocalization",
        "Total patients",
        "Total source studies",
    ]
    for disease in diseases:
        columns.extend([f"{disease} patients", f"{disease} source studies"])

    if not rows:
        return pd.DataFrame(columns=columns)

    table = pd.DataFrame(rows).sort_values(
        ["Total patients", "Colocalization"],
        ascending=[False, True],
        kind="stable",
    )
    if top_n > 0:
        table = table.head(top_n)
    table = table.reset_index(drop=True)
    table.insert(0, "Rank", np.arange(1, len(table) + 1))

    count_columns = [column for column in table if column != "Colocalization"]
    table[count_columns] = table[count_columns].astype(int)
    return table[columns]


def save_prominent_colocalizations_table(
    data: Dict[str, List[Dict]], output_path: Path, top_n: int = 10
) -> pd.DataFrame:
    """Build, save, and print the cross-disease prominent-colocalization table."""
    table = build_prominent_colocalizations_table(data, top_n=top_n)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)

    print("\nMost prominent ARG-MGE colocalizations (donor-only excluded):")
    print(table.to_string(index=False))
    print(f"Saved {output_path}")
    return table


def plot_stacked_for_disease(disease: str, records: List[Dict], mge_name_lookup: Dict[str, str]) -> None:
    pivot, top10 = build_top10_table(records)
    if pivot.empty or len(top10) == 0:
        print(f"No usable data for {disease}; skipping.")
        return

    fig, ax = plt.subplots(figsize=(11.8, 7.2), dpi=300)

    x = np.arange(len(top10))
    bottoms = np.zeros(len(top10), dtype=float)

    for status in STATUS_ORDER:
        y = pivot.loc[status].values.astype(float)
        current_bottoms = bottoms.copy()
        ax.bar(
            x,
            y,
            width=0.45,
            bottom=bottoms,
            color=STATUS_COLORS[status],
            edgecolor="black",
            linewidth=0.1,
            label=STATUS_LABELS[status],
        )

        for i, v in enumerate(y):
            if v <= 0:
                continue
            label = str(int(v)) if float(v).is_integer() else f"{v:.1f}"
            ax.text(
                x[i],
                current_bottoms[i] + (v / 2.0),
                label,
                ha="center",
                va="center",
                fontsize=14,
                fontweight="bold",
                color="white",
            )
        bottoms += y

    ax.set_xticks(x)
    ax.set_xticklabels(
        [format_colocalization_label(col, mge_name_lookup) for col in top10],
        fontsize=18,
        fontweight="bold",
        rotation=45,
        ha="right",
    )
    ax.set_ylabel("No. of patients", fontsize=20, fontweight="bold")
    # ax.set_xlabel("Most Frequent Colocalizations", fontsize=20, fontweight="bold")
    # ax.set_title(f"{disease}: Transition Types Stacked by Colocalization", fontsize=17, fontweight="bold")

    ax.tick_params(axis="y", labelsize=16)
    for t in ax.get_yticklabels():
        t.set_fontweight("bold")

    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    ax.legend(
        # title="Transition Type",
        title_fontsize=18,
        fontsize=18,
        loc="upper right",
        frameon=True,
    )

    fig.tight_layout()
    out_prefix = OUTPUT_DIR / f"{disease.lower()}_top10_colocalizations_stacked_status"
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {out_prefix.with_suffix('.png')}")
    print(f"Saved {out_prefix.with_suffix('.svg')}")


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Could not find JSON file: {DATA_PATH}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mge_name_lookup = build_mge_name_lookup(CLASSIFICATION_PATH)

    with DATA_PATH.open("r", encoding="utf-8") as f:
        data = json.load(f)

    save_prominent_colocalizations_table(data, SUMMARY_OUTPUT_PATH, top_n=10)

    for disease, records in data.items():
        plot_stacked_for_disease(disease, records, mge_name_lookup)


if __name__ == "__main__":
    main()
