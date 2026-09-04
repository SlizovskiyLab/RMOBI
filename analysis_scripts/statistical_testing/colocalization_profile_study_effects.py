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

HERE = Path(__file__).resolve().parent
ROOT_DIR = HERE.parents[1]
DATA_PATH = ROOT_DIR / "data" / "patientwise_colocalization_by_timepoint.csv"
METADATA_PATH = ROOT_DIR / "data" / "meta_data.csv"
OUTPUT_DIR = HERE / "output" / "colocalization_profile_study_effects"
FIGURE_DIR = HERE / "output" / "colocalization_profile_study_effects"
RNG_SEED = 20260803
PERMUTATIONS = 9999

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


def build_patient_colocalization_matrix(source: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    filtered = transition_rate.filter_to_snp_confirmed_args(source)
    phase = transition_rate.add_phase_presence(filtered)
    active = phase[
        phase[["Donor_p", "PreFMT_p", "PostAny_p"]].fillna(0).astype(int).any(axis=1)
    ].copy()
    active["Patient"] = pd.to_numeric(active["Patient"], errors="coerce").astype(float)
    active["Colocalization"] = (
        active["MEGARes group"].astype(str).str.strip()
        + " || "
        + active["MGE gene"].astype(str).str.strip()
    )
    active = active.drop_duplicates(["Patient", "Colocalization"])

    matrix = pd.crosstab(active["Patient"], active["Colocalization"])
    matrix = (matrix > 0).astype(int)

    patient_meta = metadata[["Patient", "disease_type", "study_data"]].drop_duplicates("Patient")
    patient_meta = patient_meta[patient_meta["Patient"].isin(matrix.index)].copy()
    patient_meta = patient_meta.sort_values(["disease_type", "study_data", "Patient"])
    matrix = matrix.reindex(patient_meta["Patient"]).fillna(0).astype(int)
    matrix.index = patient_meta["Patient"].to_numpy()
    return matrix, patient_meta


def jaccard_distance_matrix(binary_matrix: np.ndarray) -> np.ndarray:
    matrix = binary_matrix.astype(int)
    intersection = matrix @ matrix.T
    row_sums = matrix.sum(axis=1)
    union = row_sums[:, None] + row_sums[None, :] - intersection
    distances = np.where(union > 0, 1.0 - intersection / union, 0.0)
    np.fill_diagonal(distances, 0.0)
    return distances.astype(float)


def pcoa(distance_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    n = distance_matrix.shape[0]
    if n == 0:
        return np.empty((0, 0)), np.array([]), np.nan
    identity = np.eye(n)
    centering = identity - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distance_matrix ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gower)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > 1e-10
    positive_eigenvalues = eigenvalues[positive]
    coordinates = eigenvectors[:, positive] * np.sqrt(positive_eigenvalues)
    absolute_sum = np.abs(eigenvalues).sum()
    negative_fraction = float(np.abs(eigenvalues[eigenvalues < -1e-10]).sum() / absolute_sum) if absolute_sum else 0.0
    return coordinates, eigenvalues, negative_fraction


def permanova_stat(distance_matrix: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    n = len(labels)
    groups = pd.Categorical(labels)
    group_codes = groups.codes
    k = len(groups.categories)
    if n <= k or k < 2:
        return {"pseudo_f": np.nan, "r2": np.nan, "ss_between": np.nan, "ss_total": np.nan}

    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distance_matrix ** 2) @ centering
    design = np.zeros((n, k), dtype=float)
    design[np.arange(n), group_codes] = 1.0
    hat = design @ np.linalg.pinv(design.T @ design) @ design.T

    ss_total = float(np.trace(gower))
    ss_between = float(np.trace(hat @ gower))
    ss_within = ss_total - ss_between
    df_between = k - 1
    df_within = n - k
    if ss_total <= 0 or ss_within <= 0 or df_within <= 0:
        pseudo_f = np.nan
    else:
        pseudo_f = (ss_between / df_between) / (ss_within / df_within)
    r2 = ss_between / ss_total if ss_total > 0 else np.nan
    return {"pseudo_f": pseudo_f, "r2": r2, "ss_between": ss_between, "ss_total": ss_total}


def permanova(distance_matrix: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> dict[str, float]:
    observed = permanova_stat(distance_matrix, labels)
    if np.isnan(observed["pseudo_f"]):
        observed["p_value"] = np.nan
        return observed

    exceed = 0
    valid = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(labels)
        stat = permanova_stat(distance_matrix, permuted)["pseudo_f"]
        if np.isnan(stat):
            continue
        valid += 1
        exceed += int(stat >= observed["pseudo_f"] - 1e-12)
    observed["p_value"] = (exceed + 1) / (valid + 1) if valid else np.nan
    return observed


def anova_f(values: np.ndarray, labels: np.ndarray) -> float:
    groups = [values[labels == label] for label in sorted(pd.unique(labels))]
    groups = [group for group in groups if len(group)]
    if len(groups) < 2:
        return np.nan
    grand_mean = float(values.mean())
    ss_between = sum(len(group) * (float(group.mean()) - grand_mean) ** 2 for group in groups)
    ss_within = sum(float(((group - float(group.mean())) ** 2).sum()) for group in groups)
    df_between = len(groups) - 1
    df_within = len(values) - len(groups)
    if df_within <= 0 or ss_within <= 0:
        return np.nan
    return (ss_between / df_between) / (ss_within / df_within)


def dispersion_test(coordinates: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> dict[str, object]:
    if coordinates.shape[0] == 0 or coordinates.shape[1] == 0:
        return {
            "dispersion_f": np.nan,
            "dispersion_p": np.nan,
            "mean_distance_to_centroid": np.nan,
            "study_mean_distances": "",
            "note": "no positive PCoA axes available",
        }

    distances = np.zeros(coordinates.shape[0], dtype=float)
    study_means = []
    for study in sorted(pd.unique(labels)):
        mask = labels == study
        centroid = coordinates[mask].mean(axis=0)
        study_distances = np.linalg.norm(coordinates[mask] - centroid, axis=1)
        distances[mask] = study_distances
        study_means.append(f"{study}: {study_distances.mean():.4g}")

    observed = anova_f(distances, labels)
    if np.isnan(observed):
        return {
            "dispersion_f": observed,
            "dispersion_p": np.nan,
            "mean_distance_to_centroid": float(distances.mean()),
            "study_mean_distances": " | ".join(study_means),
            "note": "insufficient within-study variation for dispersion ANOVA",
        }

    exceed = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(labels)
        stat = anova_f(distances, permuted)
        exceed += int(not np.isnan(stat) and stat >= observed - 1e-12)
    return {
        "dispersion_f": observed,
        "dispersion_p": (exceed + 1) / (PERMUTATIONS + 1),
        "mean_distance_to_centroid": float(distances.mean()),
        "study_mean_distances": " | ".join(study_means),
        "note": "",
    }


def interpretation(permanova_p: float, dispersion_p: float, alpha: float = 0.05) -> str:
    permanova_sig = pd.notna(permanova_p) and permanova_p < alpha
    dispersion_sig = pd.notna(dispersion_p) and dispersion_p < alpha
    if not permanova_sig and not dispersion_sig:
        return "No detectable study separation; limited power should be acknowledged."
    if permanova_sig and not dispersion_sig:
        return "Evidence that study profile centroids differ."
    if permanova_sig and dispersion_sig:
        return "Study profile differences may reflect both centroid separation and unequal dispersion."
    return "Studies differ mainly in within-study variability."


def plot_pcoa(
    coordinates: np.ndarray,
    eigenvalues: np.ndarray,
    patient_meta: pd.DataFrame,
    disease: str,
    negative_fraction: float,
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    if coordinates.shape[1] == 0:
        return

    if coordinates.shape[1] == 1:
        coords = np.column_stack([coordinates[:, 0], np.zeros(coordinates.shape[0])])
    else:
        coords = coordinates[:, :2]

    positive = eigenvalues[eigenvalues > 1e-10]
    variance = positive / positive.sum() if positive.sum() > 0 else np.array([np.nan, np.nan])
    pc1 = float(variance[0] * 100) if len(variance) > 0 else np.nan
    pc2 = float(variance[1] * 100) if len(variance) > 1 else 0.0

    fig, ax = plt.subplots(figsize=(6.8, 5.4), dpi=300)
    studies = sorted(patient_meta["study_data"].dropna().unique())
    cmap = matplotlib.colormaps["tab20"].resampled(max(1, len(studies)))
    for idx, study in enumerate(studies):
        mask = patient_meta["study_data"].eq(study).to_numpy()
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=54,
            color=cmap(idx),
            edgecolor="white",
            linewidth=0.6,
            label=f"{study} (n={mask.sum()})",
            alpha=0.9,
        )

    ax.axhline(0, color="#dddddd", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#dddddd", linewidth=0.8, zorder=0)
    ax.set_xlabel(f"PCoA1 ({pc1:.1f}% positive-axis variance)", fontweight="bold", fontsize=12)
    ax.set_ylabel(f"PCoA2 ({pc2:.1f}% positive-axis variance)", fontweight="bold", fontsize=12)
    ax.set_title(f"{disease}", fontweight="bold", fontsize=14)
    ax.text(
        0.02,
        0.02,
        f"Negative eigenvalue fraction: {negative_fraction:.3f}",
        transform=ax.transAxes,
        fontsize=10,
        color="#444444",
    )
    ax.legend(loc="best", bbox_to_anchor=(1.02, 0.5), frameon=True, fontsize=8)
    fig.tight_layout()
    stem = disease.lower().replace(" ", "_")
    png_path = FIGURE_DIR / f"pcoa_jaccard_profile_{stem}.png"
    svg_path = FIGURE_DIR / f"pcoa_jaccard_profile_{stem}.svg"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png_path}")
    print(f"Saved {svg_path}")


def analyze_disease(matrix: pd.DataFrame, patient_meta: pd.DataFrame, disease: str, rng: np.random.Generator) -> dict[str, object]:
    disease_meta = patient_meta[patient_meta["disease_type"] == disease].copy()
    disease_matrix = matrix.loc[disease_meta["Patient"]]
    present_features = disease_matrix.columns[disease_matrix.sum(axis=0) > 0]
    disease_matrix = disease_matrix[present_features]
    labels = disease_meta["study_data"].astype(str).to_numpy()

    row = {
        "Disease": disease,
        "Patients": int(disease_matrix.shape[0]),
        "Studies": int(disease_meta["study_data"].nunique()),
        "Colocalization features": int(disease_matrix.shape[1]),
        "Mean observed features per patient": float(disease_matrix.sum(axis=1).mean()),
        "Median observed features per patient": float(disease_matrix.sum(axis=1).median()),
        "PERMANOVA pseudo-F": np.nan,
        "PERMANOVA R2": np.nan,
        "PERMANOVA P": np.nan,
        "Dispersion F": np.nan,
        "Dispersion P": np.nan,
        "Mean distance to centroid": np.nan,
        "Study mean distances to centroid": "",
        "Negative eigenvalue fraction": np.nan,
        "Interpretation": "",
        "Note": "",
    }

    if disease_matrix.shape[0] <= disease_meta["study_data"].nunique() or disease_matrix.shape[1] == 0:
        row["Note"] = "insufficient patients or nonzero colocalization features"
        return row

    distances = jaccard_distance_matrix(disease_matrix.to_numpy(dtype=int))
    coordinates, eigenvalues, negative_fraction = pcoa(distances)
    perma = permanova(distances, labels, rng)
    disper = dispersion_test(coordinates, labels, rng)
    plot_pcoa(coordinates, eigenvalues, disease_meta, disease, negative_fraction)

    row.update(
        {
            "PERMANOVA pseudo-F": perma["pseudo_f"],
            "PERMANOVA R2": perma["r2"],
            "PERMANOVA P": perma["p_value"],
            "Dispersion F": disper["dispersion_f"],
            "Dispersion P": disper["dispersion_p"],
            "Mean distance to centroid": disper["mean_distance_to_centroid"],
            "Study mean distances to centroid": disper["study_mean_distances"],
            "Negative eigenvalue fraction": negative_fraction,
            "Interpretation": interpretation(perma["p_value"], disper["dispersion_p"]),
            "Note": disper["note"],
        }
    )
    return row


def write_markdown_summary(results: pd.DataFrame, path: Path) -> None:
    lines = [
        "| Disease | Patients | Studies | Features | PERMANOVA R2 | PERMANOVA P | Dispersion P | Interpretation |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for _, row in results.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["Disease"]),
                    str(int(row["Patients"])),
                    str(int(row["Studies"])),
                    str(int(row["Colocalization features"])),
                    f"{row['PERMANOVA R2']:.4g}" if pd.notna(row["PERMANOVA R2"]) else "",
                    f"{row['PERMANOVA P']:.4g}" if pd.notna(row["PERMANOVA P"]) else "",
                    f"{row['Dispersion P']:.4g}" if pd.notna(row["Dispersion P"]) else "",
                    str(row["Interpretation"]),
                ]
            )
            + " |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    source, metadata = load_inputs()
    matrix, patient_meta = build_patient_colocalization_matrix(source, metadata)
    matrix.to_csv(OUTPUT_DIR / "patient_by_colocalization_binary_matrix.csv")
    patient_meta.to_csv(OUTPUT_DIR / "patient_profile_metadata.csv", index=False)

    rows = []
    for disease in DISEASE_ORDER:
        rows.append(analyze_disease(matrix, patient_meta, disease, rng))
    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT_DIR / "profile_study_effect_permanova_dispersion.csv", index=False)
    write_markdown_summary(results, OUTPUT_DIR / "profile_study_effect_summary.md")
    print(f"Saved profile study-effect diagnostics in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
