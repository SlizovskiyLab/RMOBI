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
from scipy.stats import kruskal


ANALYSIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = ANALYSIS_DIR.parent
DATA_PATH = ROOT_DIR / "data" / "patientwise_colocalization_by_timepoint.csv"
OUTPUT_DIR = ANALYSIS_DIR / "images" / "transition_rates"
OUTPUT_PREFIX = OUTPUT_DIR / "emergence_disappearance_log_ratio"

DISEASE_ORDER = ["MDRB", "Melanoma", "rCDI"]
POINT_COLORS = {
    "MDRB": "#8f241e",
    "Melanoma": "#3f3f3f",
    "rCDI": "#174f88",
}

MODULE_PATH = ANALYSIS_DIR / "2_transitionRate.py"
spec = importlib.util.spec_from_file_location("transition_rate", MODULE_PATH)
transition_rate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(transition_rate)
compute_patient_rates_from_definitions = transition_rate.compute_patient_rates_from_definitions


def format_p_value(value: float) -> str:
    if np.isnan(value):
        return "NA"
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def build_log_ratio_dataframe(patient_rates: pd.DataFrame) -> pd.DataFrame:
    df = patient_rates[["Disease_type", "Patient", "n_emerge", "n_disappear"]].dropna().copy()
    df = df[df["Disease_type"].isin(DISEASE_ORDER)]
    df["n_emerge"] = pd.to_numeric(df["n_emerge"], errors="coerce")
    df["n_disappear"] = pd.to_numeric(df["n_disappear"], errors="coerce")
    df = df.dropna(subset=["n_emerge", "n_disappear"])

    # Positive values should mean emergence-weighted, so emergence is the numerator.
    df["delta"] = np.log((df["n_emerge"] + 0.5) / (df["n_disappear"] + 0.5))
    return df


def epsilon_squared(h_stat: float, group_count: int, total_n: int) -> float:
    if np.isnan(h_stat) or total_n <= group_count:
        return np.nan
    return max(0.0, (h_stat - group_count + 1) / (total_n - group_count))


def run_kruskal(df: pd.DataFrame) -> tuple[float, float, float]:
    groups = [
        df.loc[df["Disease_type"] == disease, "delta"].to_numpy(dtype=float)
        for disease in DISEASE_ORDER
    ]
    if any(len(group) == 0 for group in groups):
        return np.nan, np.nan, np.nan

    result = kruskal(*groups)
    h_stat = float(result.statistic)
    effect_size = epsilon_squared(h_stat, len(groups), sum(len(group) for group in groups))
    return h_stat, float(result.pvalue), effect_size


def print_summary(df: pd.DataFrame, h_stat: float, p_value: float, effect_size: float) -> None:
    print("\nEmergence-to-disappearance log-ratio by cohort")
    print(f"Kruskal-Wallis H = {h_stat:.4f}")
    print(f"Kruskal-Wallis P = {format_p_value(p_value)}")
    print(f"Epsilon^2 effect size = {effect_size:.4f}")

    rows = []
    for disease in DISEASE_ORDER:
        values = df.loc[df["Disease_type"] == disease, "delta"].to_numpy(dtype=float)
        rows.append(
            {
                "Disease cohort": disease,
                "Patients": len(values),
                "Median delta": np.median(values) if len(values) else np.nan,
                "IQR delta": (
                    np.percentile(values, 75) - np.percentile(values, 25)
                    if len(values)
                    else np.nan
                ),
            }
        )

    summary = pd.DataFrame(rows)
    summary["Median delta"] = summary["Median delta"].map(lambda value: f"{value:.4f}")
    summary["IQR delta"] = summary["IQR delta"].map(lambda value: f"{value:.4f}")
    print(summary.to_string(index=False))


def plot_log_ratio(df: pd.DataFrame, h_stat: float, p_value: float, effect_size: float) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=300)
    positions = np.arange(1, len(DISEASE_ORDER) + 1)
    data = [
        df.loc[df["Disease_type"] == disease, "delta"].to_numpy(dtype=float)
        for disease in DISEASE_ORDER
    ]

    violins = ax.violinplot(
        data,
        positions=positions,
        widths=0.72,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )

    for body, disease in zip(violins["bodies"], DISEASE_ORDER):
        body.set_facecolor(POINT_COLORS[disease])
        body.set_edgecolor("black")
        body.set_linewidth(0.8)
        body.set_alpha(0.34)

    rng = np.random.default_rng(42)
    for pos, disease, values in zip(positions, DISEASE_ORDER, data):
        jitter = rng.normal(0, 0.055, len(values))
        ax.scatter(
            np.full(len(values), pos) + jitter,
            values,
            s=34,
            color=POINT_COLORS[disease],
            alpha=0.86,
            edgecolor="white",
            linewidth=0.45,
            zorder=3,
        )
        if len(values):
            median = float(np.median(values))
            ax.hlines(median, pos - 0.22, pos + 0.22, color="black", linewidth=2.0, zorder=4)
            ax.text(
                pos,
                ax.get_ylim()[0],
                f"n = {len(values)}",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

    ax.axhline(0, color="black", linestyle="--", linewidth=1.2, alpha=0.85)
    ax.text(
        0.02,
        0.98,
        f"Kruskal-Wallis H(2) = {h_stat:.2f}, P = {format_p_value(p_value)}\n"
        f"epsilon^2 = {effect_size:.3f}\n"
        "positive: emergence-weighted\n"
        "negative: disappearance-weighted",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="#d0d0d0", linewidth=0.8, alpha=0.9),
    )

    ax.set_xticks(positions)
    ax.set_xticklabels(DISEASE_ORDER, fontsize=12, fontweight="bold")
    # ax.set_xlabel("Cohort", fontsize=13, fontweight="bold")
    ax.set_ylabel(r"$\delta_p = \log((n_{Em}+0.5)/(n_{Dis}+0.5))$", fontsize=13, fontweight="bold")
    ax.tick_params(axis="y", labelsize=11)
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.set_axisbelow(True)

    y_min, y_max = ax.get_ylim()
    padding = (y_max - y_min) * 0.08
    ax.set_ylim(y_min - padding, y_max + padding)
    for text in ax.texts:
        if text.get_text().startswith("n = "):
            text.set_y(y_min - padding * 0.82)

    fig.tight_layout()
    fig.savefig(OUTPUT_PREFIX.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_PREFIX.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {OUTPUT_PREFIX.with_suffix('.png')}")
    print(f"Saved {OUTPUT_PREFIX.with_suffix('.svg')}")


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Could not find CSV: {DATA_PATH}")

    source = pd.read_csv(DATA_PATH)
    patient_rates = compute_patient_rates_from_definitions(source)
    ratio_df = build_log_ratio_dataframe(patient_rates)
    h_stat, p_value, effect_size = run_kruskal(ratio_df)

    print_summary(ratio_df, h_stat, p_value, effect_size)
    plot_log_ratio(ratio_df, h_stat, p_value, effect_size)


if __name__ == "__main__":
    main()
