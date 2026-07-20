#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd
os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "rmobi-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path("/tmp") / "rmobi-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata, wilcoxon


ANALYSIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = ANALYSIS_DIR.parent
DATA_PATH = ROOT_DIR / "data" / "patientwise_colocalization_by_timepoint.csv"
OUTPUT_DIR = ANALYSIS_DIR / "images" / "transition_rates"
OUTPUT_PREFIX = OUTPUT_DIR / "paired_emergence_vs_disappearance"

DISEASE_ORDER = ["MDRB", "Melanoma", "rCDI"]

MODULE_PATH = ANALYSIS_DIR / "2_transitionRate.py"
spec = importlib.util.spec_from_file_location("transition_rate", MODULE_PATH)
transition_rate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(transition_rate)
compute_patient_rates_from_definitions = transition_rate.compute_patient_rates_from_definitions


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Return BH-adjusted p-values in the original p-value order."""
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(len(p), np.nan)
    valid = np.where(~np.isnan(p))[0]
    if len(valid) == 0:
        return adjusted.tolist()

    order = valid[np.argsort(p[valid])]
    ranked_p = p[order]
    m = len(ranked_p)
    raw_adjusted = ranked_p * m / np.arange(1, m + 1)
    monotonic = np.minimum.accumulate(raw_adjusted[::-1])[::-1]
    adjusted[order] = np.clip(monotonic, 0, 1)
    return adjusted.tolist()


def paired_wilcoxon_stats(er: np.ndarray, dr: np.ndarray) -> tuple[float, float, float]:
    """Return R-style Wilcoxon V, raw p-value, and rank-biserial effect size."""
    differences = er - dr
    nonzero = differences[differences != 0]
    if len(nonzero) == 0:
        return 0.0, 1.0, 0.0

    ranks = rankdata(np.abs(nonzero))
    w_plus = float(ranks[nonzero > 0].sum())
    w_minus = float(ranks[nonzero < 0].sum())
    rank_sum = float(ranks.sum())
    effect_size = (w_plus - w_minus) / rank_sum if rank_sum else 0.0

    result = wilcoxon(er, dr, zero_method="wilcox", alternative="two-sided", mode="auto")
    return w_plus, float(result.pvalue), effect_size


def format_p_value(value: float) -> str:
    if np.isnan(value):
        return "NA"
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def summarize_cohorts(patient_rates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    raw_p_values = []

    for disease in DISEASE_ORDER:
        cohort = (
            patient_rates.loc[
                patient_rates["Disease_type"] == disease,
                ["Patient", "EmergenceRate", "DisappearanceRate"],
            ]
            .dropna()
            .copy()
        )

        er = cohort["EmergenceRate"].to_numpy(dtype=float)
        dr = cohort["DisappearanceRate"].to_numpy(dtype=float)
        v_stat, p_value, effect_size = paired_wilcoxon_stats(er, dr)
        raw_p_values.append(p_value)

        rows.append(
            {
                "Disease cohort": disease,
                "Patients": len(cohort),
                "Median ER": float(np.median(er)) if len(er) else np.nan,
                "Median DR": float(np.median(dr)) if len(dr) else np.nan,
                "Median paired difference": float(np.median(er - dr)) if len(cohort) else np.nan,
                "Wilcoxon V": v_stat,
                "Raw P": p_value,
                "Adjusted P": np.nan,
                "Effect size": effect_size,
            }
        )

    adjusted = benjamini_hochberg(raw_p_values)
    for row, adjusted_p in zip(rows, adjusted):
        row["Adjusted P"] = adjusted_p

    return pd.DataFrame(rows)


def print_summary(summary: pd.DataFrame) -> None:
    printable = summary[
        [
            "Disease cohort",
            "Patients",
            "Median ER",
            "Median DR",
            "Median paired difference",
            "Wilcoxon V",
            "Adjusted P",
            "Effect size",
        ]
    ].copy()

    for column in ["Median ER", "Median DR", "Median paired difference", "Effect size"]:
        printable[column] = printable[column].map(lambda value: f"{value:.4f}")
    printable["Wilcoxon V"] = printable["Wilcoxon V"].map(lambda value: f"{value:.2f}")
    printable["Adjusted P"] = printable["Adjusted P"].map(format_p_value)

    print("\nWithin-cohort emergence versus disappearance")
    print(printable.to_string(index=False))


def plot_paired_rates(patient_rates: pd.DataFrame, summary: pd.DataFrame) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, len(DISEASE_ORDER), figsize=(11.5, 4.0), dpi=300, sharey=True)
    colors = {"Emergence": "#b8483f", "Disappearance": "#2f6fb0"}
    rng = np.random.default_rng(42)

    for ax, disease in zip(axes, DISEASE_ORDER):
        cohort = (
            patient_rates.loc[
                patient_rates["Disease_type"] == disease,
                ["Patient", "EmergenceRate", "DisappearanceRate"],
            ]
            .dropna()
            .sort_values("Patient")
        )

        er = cohort["EmergenceRate"].to_numpy(dtype=float)
        dr = cohort["DisappearanceRate"].to_numpy(dtype=float)
        x_er = np.full(len(cohort), 0.0) + rng.normal(0, 0.025, len(cohort))
        x_dr = np.full(len(cohort), 1.0) + rng.normal(0, 0.025, len(cohort))

        for x1, x2, y1, y2 in zip(x_er, x_dr, er, dr):
            ax.plot([x1, x2], [y1, y2], color="#888888", linewidth=0.7, alpha=0.55, zorder=2)

        ax.scatter(x_er, er, s=30, color=colors["Emergence"], alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
        ax.scatter(x_dr, dr, s=30, color=colors["Disappearance"], alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)

        med_er = float(np.median(er)) if len(er) else np.nan
        med_dr = float(np.median(dr)) if len(dr) else np.nan
        ax.hlines(med_er, -0.18, 0.18, color="black", linewidth=2.0, zorder=3)
        ax.hlines(med_dr, 0.82, 1.18, color="black", linewidth=2.0, zorder=3)

        stats = summary.loc[summary["Disease cohort"] == disease].iloc[0]
        label = (
            # "Wilcoxon paired test\n"
            f"P adj = {format_p_value(float(stats['Adjusted P']))}\n"
            f"median diff = {float(stats['Median paired difference']):.3f}"
        )
        ax.text(
            0.5,
            1.28,
            label,
            ha="center",
            va="top",
            fontsize=12,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="#d0d0d0", linewidth=0.8, alpha=0.9),
        )

        ax.set_title(disease, fontsize=14, fontweight="bold")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Emergence", "Disappearance"], rotation=0, ha="center", fontsize=12, fontweight="bold")
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylim(-0.05, 1.32)
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.grid(axis="y", linestyle="--", alpha=0.25)
        ax.set_axisbelow(True)

    axes[0].set_ylabel("Patient-level transition rate", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUTPUT_PREFIX.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_PREFIX.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {OUTPUT_PREFIX.with_suffix('.png')}")
    print(f"Saved {OUTPUT_PREFIX.with_suffix('.svg')}")


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Could not find CSV: {DATA_PATH}")

    df = pd.read_csv(DATA_PATH)
    patient_rates = compute_patient_rates_from_definitions(df)
    summary = summarize_cohorts(patient_rates)
    print_summary(summary)
    plot_paired_rates(patient_rates, summary)


if __name__ == "__main__":
    main()
