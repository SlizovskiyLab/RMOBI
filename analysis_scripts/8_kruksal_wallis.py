import importlib.util
from pathlib import Path

import pandas as pd
from scipy.stats import kruskal

try:
    import scikit_posthocs as sp
except ImportError:  # pragma: no cover - optional dependency
    sp = None

MODULE_PATH = Path(__file__).resolve().with_name("2_transitionRate.py")
spec = importlib.util.spec_from_file_location("transition_rate", MODULE_PATH)
transition_rate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(transition_rate)
compute_patient_rates_from_definitions = transition_rate.compute_patient_rates_from_definitions

# Metrics to test
METRICS = [
    "PersistenceRate",
    "DisappearanceRate",
    "EmergenceRate",
    "TransferShare",
]

DISEASE_ORDER = ["MDRB", "Melanoma", "rCDI"]
DATA_PATH = Path(__file__).resolve().parents[1] / "data" / "patientwise_colocalization_by_timepoint.csv"

def run_kruskal_tests(patient_rates):
    """
    Performs Kruskal-Wallis tests across disease cohorts for each
    transition-rate metric. If significant, performs Dunn's post-hoc test.
    """

    for metric in METRICS:

        print("=" * 70)
        print(metric)

        groups = []

        for disease in DISEASE_ORDER:
            vals = (
                patient_rates.loc[
                    patient_rates["Disease_type"] == disease,
                    metric
                ]
                .dropna()
                .values
            )
            groups.append(vals)

        # Overall Kruskal-Wallis
        H, p = kruskal(*groups)
        n = sum(len(g) for g in groups)
        k = len(groups)

        eps2 = epsilon_squared(H, k, n)

        print(f"Kruskal-Wallis H = {H:.4f}")
        print(f"p-value            = {p:.6f}")
        print(f"Epsilon^2 = {eps2:.3f}")
        print(f"Effect size = {interpret_effect_size(eps2)}")

        if p < 0.05:
            if sp is None:
                print("\nscikit-posthocs is not installed; skipping Dunn post-hoc test.")
            else:
                print("\nPost-hoc Dunn test (Holm correction):")

                dunn = sp.posthoc_dunn(
                    patient_rates,
                    val_col=metric,
                    group_col="Disease_type",
                    p_adjust="holm",
                )

                dunn = dunn.loc[DISEASE_ORDER, DISEASE_ORDER]

                print(dunn.round(4))

        else:
            print("\nNo significant overall difference.")



def epsilon_squared(H, k, n):
    """
    Epsilon-squared effect size for Kruskal-Wallis.

    Parameters
    ----------
    H : float
        Kruskal-Wallis H statistic.
    k : int
        Number of groups.
    n : int
        Total number of observations.

    Returns
    -------
    float
    """
    return max(0, (H - k + 1) / (n - k))


def interpret_effect_size(eps):
    if eps < 0.01:
        return "negligible"
    elif eps < 0.08:
        return "small"
    elif eps < 0.26:
        return "moderate"
    else:
        return "large"


if __name__ == "__main__":

    df = pd.read_csv(DATA_PATH)

    patient_rates = compute_patient_rates_from_definitions(df)

    run_kruskal_tests(patient_rates)