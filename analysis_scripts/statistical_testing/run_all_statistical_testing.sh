#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v python3 >/dev/null || { echo "python3 is required" >&2; exit 1; }
command -v Rscript >/dev/null || { echo "Rscript is required" >&2; exit 1; }

echo "Python: prepare model inputs"
python3 "$SCRIPT_DIR/prepare_emergence_disappearance_beta_binomial.py"

echo "Python: colocalization profile study effects"
python3 "$SCRIPT_DIR/colocalization_profile_study_effects.py"

echo "Python: random-effects meta-analysis"
python3 "$SCRIPT_DIR/random_effects_meta_analysis.py"

echo "R: Model 1 (all ARGs)"
Rscript "$SCRIPT_DIR/fit_emergence_disappearance_beta_binomial.R"

# echo "R: Model 1 (clinically relevant ARGs)"
# Rscript "$SCRIPT_DIR/fit_emergence_disappearance_beta_binomial.R" \
#   "$SCRIPT_DIR/output/model_1_clinically_relevant/emergence_disappearance_patients.csv" \
#   "$SCRIPT_DIR/output/model_1_clinically_relevant"

echo "R: transition-type models"
Rscript "$SCRIPT_DIR/fit_transition_type_beta_binomial.R"

echo "R: MGE-class interaction model"
Rscript "$SCRIPT_DIR/fit_mge_class_interaction_beta_binomial.R"

echo "R: leave-one-study-out analysis"
Rscript "$SCRIPT_DIR/fit_leave_one_study_out.R"

echo "R: Model 1 figures"
Rscript "$SCRIPT_DIR/plot_emergence_disappearance_balance.R"

echo "All statistical-testing analyses completed."
