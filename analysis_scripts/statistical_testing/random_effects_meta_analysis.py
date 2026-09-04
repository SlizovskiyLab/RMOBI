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
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, norm, t


HERE = Path(__file__).resolve().parent
ROOT_DIR = HERE.parents[1]
DATA_PATH = ROOT_DIR / "data" / "patientwise_colocalization_by_timepoint.csv"
METADATA_PATH = ROOT_DIR / "data" / "meta_data.csv"
OUTPUT_DIR = HERE / "output" / "random_effects_meta_analysis"
FIGURE_DIR = HERE / "output" / "random_effects_meta_analysis"

DISEASE_ORDER = ["MDRB", "Melanoma", "rCDI"]
DISEASE_COLORS = {
    "MDRB": "#8f241e",
    "Melanoma": "#4a4a4a",
    "rCDI": "#1f5f99",
}

MODULE_PATH = ROOT_DIR / "analysis_scripts" / "2_transitionRate.py"
spec = importlib.util.spec_from_file_location("transition_rate", MODULE_PATH)
transition_rate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(transition_rate)
compute_patient_rates_from_definitions = transition_rate.compute_patient_rates_from_definitions


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Could not find main CSV: {DATA_PATH}")
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"Could not find patient metadata CSV: {METADATA_PATH}")

    source = pd.read_csv(DATA_PATH)
    metadata = pd.read_csv(METADATA_PATH)
    metadata["Patient"] = pd.to_numeric(metadata["Patient"], errors="coerce").astype("Int64")
    metadata = metadata.dropna(subset=["Patient", "disease_type", "study_data"]).copy()
    metadata["Patient"] = metadata["Patient"].astype(float)
    return source, metadata


def build_patient_log_ratios(source: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    patient_rates = compute_patient_rates_from_definitions(source)
    patient_rates["Patient"] = pd.to_numeric(patient_rates["Patient"], errors="coerce").astype(float)
    patient_study = metadata[["Patient", "study_data"]].drop_duplicates("Patient")
    df = patient_rates.merge(patient_study, on="Patient", how="left")
    df = df.dropna(subset=["study_data", "n_emerge", "n_disappear"]).copy()
    df = df[df["Disease_type"].isin(DISEASE_ORDER)].copy()
    df["n_emerge"] = pd.to_numeric(df["n_emerge"], errors="coerce")
    df["n_disappear"] = pd.to_numeric(df["n_disappear"], errors="coerce")
    df["EmergenceDisappearanceLogRatio"] = np.log(
        (df["n_emerge"] + 0.5) / (df["n_disappear"] + 0.5)
    )
    return df.sort_values(["Disease_type", "study_data", "Patient"])


def study_specific_estimates(patient_ratios: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (disease, study), group in patient_ratios.groupby(["Disease_type", "study_data"], observed=True):
        values = group["EmergenceDisappearanceLogRatio"].dropna().to_numpy(dtype=float)
        if len(values) == 0:
            continue
        sd = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
        se = sd / np.sqrt(len(values)) if len(values) > 1 else np.nan
        rows.append(
            {
                "Disease": disease,
                "Study": study,
                "Patients": int(len(values)),
                "Estimate": float(np.mean(values)),
                "SD": sd,
                "SE": se,
                "Variance": se**2 if pd.notna(se) else np.nan,
                "Median": float(np.median(values)),
                "IQR": float(np.percentile(values, 75) - np.percentile(values, 25)),
                "Minimum": float(np.min(values)),
                "Maximum": float(np.max(values)),
                "Emergence events": int(group["n_emerge"].sum()),
                "Disappearance events": int(group["n_disappear"].sum()),
            }
        )
    estimates = pd.DataFrame(rows)
    estimates["Disease"] = pd.Categorical(estimates["Disease"], categories=DISEASE_ORDER, ordered=True)
    return estimates.sort_values(["Disease", "Study"]).reset_index(drop=True)


def fixed_effect_heterogeneity(yi: np.ndarray, vi: np.ndarray) -> dict[str, float]:
    weights = 1.0 / vi
    fixed_mu = float(np.sum(weights * yi) / np.sum(weights))
    q_stat = float(np.sum(weights * (yi - fixed_mu) ** 2))
    df = len(yi) - 1
    p_q = float(chi2.sf(q_stat, df)) if df > 0 else np.nan
    i2 = max(0.0, (q_stat - df) / q_stat) * 100.0 if q_stat > 0 and df > 0 else 0.0
    return {"fixed_mu": fixed_mu, "Q": q_stat, "Q_df": df, "P_Q": p_q, "I2_percent": i2}


def reml_objective(tau2: float, yi: np.ndarray, vi: np.ndarray) -> float:
    total_variance = vi + tau2
    weights = 1.0 / total_variance
    mu = np.sum(weights * yi) / np.sum(weights)
    residual = yi - mu
    return float(
        np.sum(np.log(total_variance))
        + np.log(np.sum(weights))
        + np.sum(weights * residual**2)
    )


def estimate_tau2_reml(yi: np.ndarray, vi: np.ndarray) -> float:
    heterogeneity = fixed_effect_heterogeneity(yi, vi)
    q_stat = heterogeneity["Q"]
    df = heterogeneity["Q_df"]
    weights = 1.0 / vi
    denominator = np.sum(weights) - np.sum(weights**2) / np.sum(weights)
    dl_tau2 = max(0.0, (q_stat - df) / denominator) if denominator > 0 else 0.0
    upper = max(float(np.var(yi, ddof=1)) * 10.0, dl_tau2 * 10.0, 1.0)

    result = minimize_scalar(
        reml_objective,
        args=(yi, vi),
        bounds=(0.0, upper),
        method="bounded",
        options={"xatol": 1e-10},
    )
    if not result.success:
        return dl_tau2
    return max(0.0, float(result.x))


def random_effects_meta_analysis(estimates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for disease in DISEASE_ORDER:
        disease_estimates = estimates[estimates["Disease"].astype(str) == disease].dropna(subset=["Variance"]).copy()
        yi = disease_estimates["Estimate"].to_numpy(dtype=float)
        vi = disease_estimates["Variance"].to_numpy(dtype=float)
        valid = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
        yi = yi[valid]
        vi = vi[valid]
        k = len(yi)

        row = {
            "Disease": disease,
            "Studies": int(k),
            "Patients": int(disease_estimates.loc[valid, "Patients"].sum()) if k else 0,
            "Pooled estimate": np.nan,
            "CI lower": np.nan,
            "CI upper": np.nan,
            "SE pooled": np.nan,
            "HK scale": np.nan,
            "tau2": np.nan,
            "I2 percent": np.nan,
            "Q": np.nan,
            "Q df": np.nan,
            "P_Q": np.nan,
            "Method": "REML random-effects with Hartung-Knapp CI",
            "Note": "",
        }
        if k < 2:
            row["Note"] = "fewer than two studies with estimable within-study SE"
            rows.append(row)
            continue

        tau2 = estimate_tau2_reml(yi, vi)
        weights = 1.0 / (vi + tau2)
        pooled = float(np.sum(weights * yi) / np.sum(weights))
        standard_se = float(np.sqrt(1.0 / np.sum(weights)))
        hk_q = float(np.sum(weights * (yi - pooled) ** 2) / (k - 1))
        hk_scale = max(1.0, hk_q)
        hk_se = standard_se * np.sqrt(hk_scale)
        critical = float(t.ppf(0.975, k - 1))
        heterogeneity = fixed_effect_heterogeneity(yi, vi)

        row.update(
            {
                "Pooled estimate": pooled,
                "CI lower": pooled - critical * hk_se,
                "CI upper": pooled + critical * hk_se,
                "SE pooled": hk_se,
                "HK scale": hk_scale,
                "tau2": tau2,
                "I2 percent": heterogeneity["I2_percent"],
                "Q": heterogeneity["Q"],
                "Q df": heterogeneity["Q_df"],
                "P_Q": heterogeneity["P_Q"],
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def add_study_confidence_intervals(estimates: pd.DataFrame) -> pd.DataFrame:
    estimates = estimates.copy()
    estimates["CI lower"] = estimates["Estimate"] - norm.ppf(0.975) * estimates["SE"]
    estimates["CI upper"] = estimates["Estimate"] + norm.ppf(0.975) * estimates["SE"]
    return estimates


def forest_xlim(estimates: pd.DataFrame, pooled: pd.DataFrame) -> tuple[float, float]:
    bounds = pd.concat(
        [
            estimates[["CI lower", "CI upper"]].stack(),
            pooled[["CI lower", "CI upper"]].stack(),
        ]
    ).dropna()
    if bounds.empty:
        return (-2.0, 2.0)
    lo = float(bounds.min())
    hi = float(bounds.max())
    pad = max(0.4, (hi - lo) * 0.16)
    return lo - pad, hi + pad


def plot_forest(estimates: pd.DataFrame, pooled: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(DISEASE_ORDER),
        1,
        figsize=(9.0, 9.8),
        dpi=300,
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 0.9, 1.8]},
    )
    x_min, x_max = forest_xlim(estimates, pooled)

    for ax, disease in zip(axes, DISEASE_ORDER):
        disease_estimates = estimates[estimates["Disease"].astype(str) == disease].copy()
        disease_pooled = pooled[pooled["Disease"] == disease].iloc[0]
        y_positions = np.arange(len(disease_estimates), 0, -1)
        color = DISEASE_COLORS[disease]

        for y, (_, row) in zip(y_positions, disease_estimates.iterrows()):
            ci_low = row["CI lower"]
            ci_high = row["CI upper"]
            estimate = row["Estimate"]
            if pd.notna(ci_low) and pd.notna(ci_high):
                ax.hlines(y, ci_low, ci_high, color="#333333", linewidth=1.2)
                ax.plot([ci_low, ci_high], [y, y], marker="|", color="#333333", linestyle="None", markersize=8)
            ax.scatter(estimate, y, s=46 + 8 * row["Patients"], color=color, edgecolor="white", linewidth=0.7, zorder=3)
            right_text = (
                f"{estimate:.2f} [{ci_low:.2f}, {ci_high:.2f}]"
                if pd.notna(ci_low) and pd.notna(ci_high)
                else f"{estimate:.2f} [SE unavailable]"
            )
            ax.text(x_max, y, right_text, ha="right", va="center", fontsize=10)

        pooled_y = 0
        if pd.notna(disease_pooled["Pooled estimate"]):
            center = disease_pooled["Pooled estimate"]
            ci_low = disease_pooled["CI lower"]
            ci_high = disease_pooled["CI upper"]
            diamond_y = np.array([pooled_y, pooled_y + 0.22, pooled_y, pooled_y - 0.22])
            diamond_x = np.array([ci_low, center, ci_high, center])
            ax.fill(diamond_x, diamond_y, color=color, alpha=0.82, edgecolor="black", linewidth=0.8)
            ax.text(
                x_max,
                pooled_y,
                f"{center:.2f} [{ci_low:.2f}, {ci_high:.2f}]",
                ha="right",
                va="center",
                fontsize=10,
                fontweight="bold",
            )

        labels = [f"{row['Study']} (n={int(row['Patients'])})" for _, row in disease_estimates.iterrows()]
        labels.append("Pooled")
        ax.set_yticks(list(y_positions) + [pooled_y])
        ax.set_yticklabels(labels, fontsize=10)
        ax.axvline(0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
        ax.set_xlim(x_min, x_max)
        ax.set_title(
            f"{disease}",
            loc="center",
            fontsize=11,
            fontweight="bold",
        )
        ax.text(
            0.01,
            0.02,
            f"I2={disease_pooled['I2 percent']:.1f}%, tau2={disease_pooled['tau2']:.3g}, "
            f"Q={disease_pooled['Q']:.2f}, P_Q={disease_pooled['P_Q']:.3g}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=10,
            color="#333333",
        )
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)

    axes[-1].set_xlabel("Emergence-to-disappearance log-ratio\npositive = emergence-weighted")
    fig.suptitle("Study-specific and pooled emergence/disappearance balance", fontsize=14, fontweight="bold")
    fig.tight_layout()
    png_path = FIGURE_DIR / "random_effects_emergence_disappearance_forest.png"
    svg_path = FIGURE_DIR / "random_effects_emergence_disappearance_forest.svg"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png_path}")
    print(f"Saved {svg_path}")


def write_markdown_summary(pooled: pd.DataFrame, path: Path) -> None:
    lines = [
        "| Disease | Studies | Patients | Pooled estimate | 95% CI | I2 (%) | tau2 | Q | P_Q |",
        "| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for _, row in pooled.iterrows():
        ci = (
            f"[{row['CI lower']:.4g}, {row['CI upper']:.4g}]"
            if pd.notna(row["CI lower"]) and pd.notna(row["CI upper"])
            else ""
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["Disease"]),
                    str(int(row["Studies"])),
                    str(int(row["Patients"])),
                    f"{row['Pooled estimate']:.4g}" if pd.notna(row["Pooled estimate"]) else "",
                    ci,
                    f"{row['I2 percent']:.4g}" if pd.notna(row["I2 percent"]) else "",
                    f"{row['tau2']:.4g}" if pd.notna(row["tau2"]) else "",
                    f"{row['Q']:.4g}" if pd.notna(row["Q"]) else "",
                    f"{row['P_Q']:.4g}" if pd.notna(row["P_Q"]) else "",
                ]
            )
            + " |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    source, metadata = load_inputs()
    patient_ratios = build_patient_log_ratios(source, metadata)
    study_estimates = add_study_confidence_intervals(study_specific_estimates(patient_ratios))
    pooled_estimates = random_effects_meta_analysis(study_estimates)

    patient_ratios.to_csv(OUTPUT_DIR / "patient_emergence_disappearance_log_ratios.csv", index=False)
    study_estimates.to_csv(OUTPUT_DIR / "study_specific_log_ratio_estimates.csv", index=False)
    pooled_estimates.to_csv(OUTPUT_DIR / "random_effects_meta_analysis_summary.csv", index=False)
    write_markdown_summary(pooled_estimates, OUTPUT_DIR / "random_effects_meta_analysis_summary.md")
    plot_forest(study_estimates, pooled_estimates)
    print(f"Saved random-effects meta-analysis outputs in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
