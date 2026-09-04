# Emergence/disappearance beta-binomial mixed model

Python constructs and checks the CSVs; R fits the model. Neither script invokes the other.

## Run all statistical-testing and cross study assessment analysis

Use the shell runner from any working directory:

```bash
bash analysis_scripts/statistical_testing/run_all_statistical_testing.sh
```

The runner keeps the languages separate and stops immediately if any command fails. It runs all
Python work first, in below order:

1. prepare the CSV inputs for the GLMMs and LOSO analysis;
2. analyze colocalization-profile study effects;
3. run the Python random-effects meta-analysis.

It then runs the R analyses in below order:

1. Model 1 for all ARGs;
2. the persistence, disappearance, emergence, and transfer models;
3. the cohort-by-MGE-class interaction model;
4. the leave-one-study-out analysis;
5. the Model 1 figures.

This shell script is preferable to separate Python and R runners because there is one entry point
and one enforced dependency order, while Python still never invokes R. The individual commands
below remain available when only one analysis needs to be rerun.

## 1. Run Python

From the repository root:

```bash
python3 analysis_scripts/statistical_testing/prepare_emergence_disappearance_beta_binomial.py
```

The patient CSV contains one row per patient represented after the established SNP filter:
`Patient`, `cohort`, `study`, `n_Em`, `n_Dis`, `k_post`, and median recipient values for
`num_contigs`, `total_length`, and `non_host_length`. An event-level audit CSV is also exported.

Patient cohort, study, and `k_post` all come from `data/meta_data.csv`. Cohort and study use
`disease_type` and `study_data`; `k_post` is the number of unique sample IDs whose normalized
`donor_pre_post` value is `PostFMT`. `patient_metadata.csv` is not used.

Sample-level `num_contigs`, `total_length`, and `non_host_length` come from `data/meta_data.csv`.
Python takes the median of each metric across all recipient samples for each patient. R applies
`log1p` and standardizes each metric before model fitting.

## 2. Install R packages once

In an R console:

```r
install.packages(c("glmmTMB", "emmeans", "DHARMa"))
```

## 3. Run R separately

From the repository root:

```bash
Rscript "analysis_scripts/statistical_testing/fit_emergence_disappearance_beta_binomial.R"
```

# Model 2: Emergence-to-dissappearance balance
The default fitted model is:

```r
glmmTMB(cbind(n_Em, n_Dis) ~ cohort + log(k_post) +
          log_num_contigs_z + log_total_length_z + log_non_host_length_z +
          (1 | study),
        ziformula = ~0, family = betabinomial(link = "logit"),
        data = pat, REML = FALSE)
```
`log(k_post)` and both patient-level assembly metrics are always included in Model 1.

Model 1 has one aggregated row per patient, so `(1 | Patient)` is not identifiable separately
from beta-binomial dispersion. In a diagnostic fit its SD collapsed to approximately zero and
its addition worsened both AIC and BIC; it is therefore not included. A patient random intercept
is appropriate in Model 4, where each patient contributes multiple MGE-class rows.

A zero-inflated beta-binomial candidate (`ziformula = ~1`) was also checked. It had the same
log-likelihood as the non-zero-inflated model, an effectively absent zero-inflation component,
and worse AIC/BIC. The primary model therefore explicitly uses `ziformula = ~0`. This cannot
make rows with `n_Em + n_Dis == 0` valid binomial observations; those rows remain excluded.

The script tests the primary model with `DHARMa::testZeroInflation(..., alternative = "greater")`.
Only if this test detects excess zero-emergence outcomes at `p < 0.05` does it fit a secondary
`ziformula = ~1` candidate and export its AIC/BIC comparison. The primary fitted object and
reported cohort contrasts always come from `ziformula = ~0`.

The full, reduced, and parsimonious models all use maximum likelihood (`REML = FALSE`), allowing
the nested likelihood-ratio, AIC, and BIC comparisons to use the fitted model objects directly.


The script also fits a nested reduced sensitivity model excluding `log_num_contigs_z`,
`log_total_length_z`, and `log_non_host_length_z`. It exports `anova(reduced_fit, fit)` as their joint likelihood-ratio test and a table of both models' log-likelihood, AIC, BIC, and convergence diagnostics. The full
assembly-adjusted model remains primary; individual covariate p-values alone are not used for
post-hoc confounder selection.

An intermediate parsimonious sensitivity model retains `log_non_host_length_z` but excludes
`log_num_contigs_z` and `log_total_length_z`. Its serialized fit, complete summary, and its
ML likelihood-ratio/AIC/BIC comparison with the full model are exported separately.

R excludes and lists patients with `n_Em + n_Dis == 0`, because a binomial response has no
trials for those rows. It reports Holm-adjusted odds ratios and 95% confidence intervals for
MDRB versus rCDI and melanoma versus rCDI, the study random-intercept SD, convergence
diagnostics, the model summary, and a serialized model.


# Model 2: PR, DR, ER and transfer

From the repository root, create the four minimal model inputs with the existing Python script:

```bash
python3 analysis_scripts/statistical_testing/prepare_emergence_disappearance_beta_binomial.py
```

Then run Model 2 separately in R:

```bash
Rscript analysis_scripts/statistical_testing/fit_transition_type_beta_binomial.R
```

Each model uses `n_total = persistence + disappearance + emergence + transfer`. The transfer
CSV is restricted to patients with a donor sample in `data/meta_data.csv`. Holm correction is
applied across the four transition models separately for each cohort contrast.

For each transition, R writes a complete `<transition>_model_summary.txt` with fixed effects,
beta-binomial dispersion, random-effect variance/SD, convergence information, and model-fit
statistics. `transition_model_fit_statistics.csv` consolidates sample size, log-likelihood,
AIC, BIC, −2 log-likelihood, residual degrees of freedom, dispersion, and the approximate Pearson
chi-square/df goodness-of-fit diagnostic.


# Model 3: cohort × MGE-class interaction

```bash
python3 analysis_scripts/statistical_testing/prepare_emergence_disappearance_beta_binomial.py
Rscript analysis_scripts/statistical_testing/fit_mge_class_interaction_beta_binomial.R
```

The Model 3 CSV includes the patient-level median recipient `num_contigs`, `total_length`, and
`non_host_length` values. Three nested interaction models are fitted by maximum likelihood:

- reduced: `cohort * MGE_class + log(k_post)`;
- parsimonious: reduced plus standardized log-transformed `non_host_length`;
- full: reduced plus all three standardized log-transformed assembly covariates.

All three include `(1|study) + (1|patient)`. Their sequential likelihood-ratio, AIC, and BIC
comparison is exported to `assembly_covariate_model_comparison.csv`. A separate additive null
model containing all three assembly covariates is compared with the full interaction model, so
the primary test remains specifically `cohort × MGE_class = 0`. MGE classes absent from one or
more cohorts are excluded because their cohort interactions cannot be estimated. ICE rows are
excluded before all Model 3 fits and written to `excluded_ice_rows.csv` for audit. Results are
saved under `output/model_3`.
Post-hoc `emmeans` comparisons test MDRB vs rCDI and melanoma vs rCDI within each of the four
retained non-ICE MGE classes; their p-values are Holm-corrected together across all eight comparisons.
The adjusted emergence-probability plot shows response-scale `emmeans` and 95% CIs for rCDI,
MDRB, and melanoma within each MGE class, omitting sparse ICE estimates.


# Leave-one-study-out sensitivity analysis

```bash
python3 analysis_scripts/statistical_testing/leave_one_study_out.py
Rscript analysis_scripts/statistical_testing/fit_leave_one_study_out.R
```

The full model and every leave-one-study-out iteration use
`cbind(n_Em, n_Dis) ~ cohort + log(k_post) + (1 | study)` without sequencing depth. Outputs
include cohort ORs, 95% CIs, p-values, study SD, dispersion, convergence diagnostics, and a
single-panel LOSO coefficient forest plot.