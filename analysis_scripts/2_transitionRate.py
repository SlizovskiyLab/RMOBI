import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


DISEASE_ORDER = ["MDRB", "Melanoma", "rCDI"]
DATA_PATH = Path("../data/patientwise_colocalization_by_timepoint.csv")
OUTPUT_DIR = Path("images/transition_rates/")
# ---------------------------
#
# ---------------------------
def postfmt_bin(colname: str):
    m = re.fullmatch(r"PostFMT_(\d{3})", colname)
    if not m:
        return None
    day = int(m.group(1))
    if 1 <= day <= 30:
        return "PostFMT1"
    if 31 <= day <= 60:
        return "PostFMT2"
    return "PostFMT3"

def add_phase_presence(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    post_cols = [c for c in df.columns if postfmt_bin(c) is not None]
    groups = {"PostFMT1": [], "PostFMT2": [], "PostFMT3": []}
    for c in post_cols:
        groups[postfmt_bin(c)].append(c)

    df["Donor_p"]  = (pd.to_numeric(df["Donor"],  errors="coerce").fillna(0) > 0).astype(int)
    df["PreFMT_p"] = (pd.to_numeric(df["PreFMT"], errors="coerce").fillna(0) > 0).astype(int)

    for g in ["PostFMT1", "PostFMT2", "PostFMT3"]:
        if groups[g]:
            block = df[groups[g]].apply(pd.to_numeric, errors="coerce").fillna(0)
            df[f"{g}_p"] = (block.max(axis=1) > 0).astype(int)
        else:
            df[f"{g}_p"] = 0

    df["PostAny_p"] = (df["PostFMT1_p"] | df["PostFMT2_p"] | df["PostFMT3_p"]).astype(int)
    return df

def compute_patient_rates_from_definitions(
    df: pd.DataFrame,
    coloc_id_cols=("Patient", "MEGARes group", "MGE gene"),
):
    df = add_phase_presence(df)
    df_u = df.drop_duplicates(list(coloc_id_cols)).copy()

    D = df_u["Donor_p"].astype(bool)
    P = df_u["PreFMT_p"].astype(bool)
    T = df_u["PostAny_p"].astype(bool)

    df_u["Persistence"]    = (P & T).astype(int)
    df_u["Disappearance"]  = (P & ~T).astype(int)
    df_u["Transfer"]       = (D & ~P & T).astype(int)
    df_u["Emergence"]      = ((~D) & (~P) & T).astype(int)

    df_u["RecipientObserved"] = (P | T).astype(int)

    g = df_u.groupby(["Disease_type", "Patient"], as_index=False).agg(
        n_recipient_obs=("RecipientObserved", "sum"),
        n_persist=("Persistence", "sum"),
        n_disappear=("Disappearance", "sum"),
        n_emerge=("Emergence", "sum"),
        n_transfer=("Transfer", "sum"),
    )

    def safe_div(a, b):
        return np.where(b > 0, a / b, np.nan)

    g["PersistenceRate"]    = safe_div(g["n_persist"],   g["n_recipient_obs"])
    g["DisappearanceRate"]  = safe_div(g["n_disappear"], g["n_recipient_obs"])
    g["EmergenceRate"]      = safe_div(g["n_emerge"],    g["n_recipient_obs"])
    g["TransferShare"]      = safe_div(g["n_transfer"],  g["n_recipient_obs"])

    return g


# ---------------------------
# Plotting with colors + naming
# ---------------------------
METRICS = ["PersistenceRate", "DisappearanceRate", "EmergenceRate", "TransferShare"]
METRIC_LABELS = {
    "PersistenceRate": "PR ",
    "DisappearanceRate": "DR ",
    "EmergenceRate": "ER ",
    "TransferShare": "TR ",   # rename as requested
}

# Red color scale for boxplots
RED_COLORS = plt.cm.Reds(np.linspace(0.45, 0.9, len(METRICS)))

def plot_one_disease_4_boxplots(patient_rates, disease, metrics=METRICS,
                                figsize=(6.2, 4.8), dpi=240,
                                show_points=True):
    d = patient_rates[patient_rates["Disease_type"] == disease].copy()
    if d.empty:
        print(f"No data for {disease}; skipping.")
        return

    data = [d[m].dropna().values for m in metrics]
    labels = [METRIC_LABELS.get(m, m) for m in metrics]

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    bp = ax.boxplot(
        data,
        labels=labels,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.6),
        boxprops=dict(edgecolor="black", linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=1.2),
        capprops=dict(color="black", linewidth=1.2),
    )

    for b in bp["boxes"]:
        b.set_facecolor("white")

    for b, c in zip(bp["boxes"], RED_COLORS):
        b.set_facecolor(c)

    # Color each set of points to match its metric
    if show_points:
        for i, (vals, c) in enumerate(zip(data, RED_COLORS), start=1):
            if len(vals) == 0:
                continue
            x = np.random.normal(loc=i, scale=0.05, size=len(vals))
            ax.scatter(
                x, vals,
                s=16, alpha=0.6,
                edgecolors="none",
                c=[c]
            )

    ax.set_ylabel("Rate", fontsize=20, fontweight="bold")

    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_fontweight("bold")

    ax.set_ylim(0, 1)
    fig.tight_layout()
    out_file = OUTPUT_DIR / f"boxplot_transition_rates_{disease.lower()}.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    fig.savefig(out_file.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_file}")


def plot_all_diseases_separately(patient_rates, metrics=METRICS):
    diseases = [d for d in DISEASE_ORDER if d in patient_rates["Disease_type"].dropna().unique()]
    for dis in diseases:
        plot_one_disease_4_boxplots(patient_rates, dis, metrics=metrics)


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Could not find CSV: {DATA_PATH}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA_PATH)
    patient_rates = compute_patient_rates_from_definitions(df)
    plot_all_diseases_separately(patient_rates)


if __name__ == "__main__":
    main()
