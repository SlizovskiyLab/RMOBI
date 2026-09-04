#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
})
if (!requireNamespace("DHARMa", quietly = TRUE)) {
  stop("Package 'DHARMa' is required. Install it once with install.packages('DHARMa').")
}

args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_flag)) dirname(normalizePath(sub("^--file=", "", script_flag[1]))) else getwd()
default_input <- file.path(script_dir, "output", "model_1", "emergence_disappearance_patients.csv")
input_csv <- if (length(args) >= 1) args[1] else default_input
output_dir <- if (length(args) >= 2) args[2] else dirname(input_csv)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


pat_all <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
assembly_columns <- c(
  "median_recipient_num_contigs", "median_recipient_total_length",
  "median_recipient_non_host_length"
)
required <- c("Patient", "cohort", "study", "n_Em", "n_Dis", "k_post", assembly_columns)
missing_columns <- setdiff(required, names(pat_all))
if (length(missing_columns)) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))

# Convert columns to appropriate types and check for missing or invalid values.
pat_all$cohort <- factor(pat_all$cohort, levels = c("rCDI", "MDRB", "Melanoma"))
pat_all$study <- factor(pat_all$study)
pat_all$n_Em <- as.integer(pat_all$n_Em)
pat_all$n_Dis <- as.integer(pat_all$n_Dis)
pat_all$k_post <- as.integer(pat_all$k_post)
if (anyNA(pat_all[, required])) stop("Required model columns contain missing values.")
if (any(pat_all$n_Em < 0 | pat_all$n_Dis < 0)) stop("Event counts cannot be negative.")
if (any(pat_all$k_post <= 0)) stop("k_post must be positive before log(k_post) is used.")

# Assembly metrics are patient medians across recipient samples. Log-transform
# their skewed distributions and standardize so coefficients are per 1 SD.
for (column in assembly_columns) {
  pat_all[[column]] <- as.numeric(pat_all[[column]])
  if (anyNA(pat_all[[column]]) || any(pat_all[[column]] <= 0)) {
    stop(column, " must contain positive, nonmissing values.")
  }
}
pat_all$log_num_contigs_z <- as.numeric(scale(log1p(pat_all$median_recipient_num_contigs)))
pat_all$log_total_length_z <- as.numeric(scale(log1p(pat_all$median_recipient_total_length)))
pat_all$log_non_host_length_z <- as.numeric(scale(log1p(pat_all$median_recipient_non_host_length)))

# Subset the data to include only patients with at least one trial (n_Em + n_Dis > 0) and save any excluded patients to a CSV file.
pat <- subset(pat_all, n_Em + n_Dis > 0)
excluded <- subset(pat_all, n_Em + n_Dis == 0)
if (nrow(excluded)) {
  write.csv(excluded, file.path(output_dir, "patients_excluded_zero_trials.csv"), row.names = FALSE)
  warning(nrow(excluded), " patient(s) excluded because n_Em + n_Dis = 0.")
}

# Fit Model 1 with mandatory sampling-frequency and assembly covariates.
model_formula <- cbind(n_Em, n_Dis) ~ cohort + log(k_post) +
  log_num_contigs_z + log_total_length_z + log_non_host_length_z + (1 | study)
fit <- glmmTMB(
  model_formula,
  ziformula = ~0,
  family = betabinomial(link = "logit"),
  data = pat,
  REML = FALSE
)
saveRDS(fit, file.path(output_dir, "beta_binomial_model.rds"))
capture.output(summary(fit), file = file.path(output_dir, "beta_binomial_model_summary.txt"))

# Compare the full assembly-adjusted model with its nested reduced model.
reduced_fit <- glmmTMB(
  cbind(n_Em, n_Dis) ~ cohort + log(k_post) + (1 | study),
  ziformula = ~0,
  family = betabinomial(link = "logit"),
  data = pat,
  REML = FALSE
)
saveRDS(reduced_fit, file.path(output_dir, "beta_binomial_model_reduced.rds"))
capture.output(summary(reduced_fit),
               file = file.path(output_dir, "beta_binomial_model_reduced_summary.txt"))
assembly_lrt <- anova(reduced_fit, fit)
capture.output(assembly_lrt,
               file = file.path(output_dir, "assembly_covariate_likelihood_ratio_test.txt"))
fit_comparison <- data.frame(
  model = c("Reduced: cohort + log(k_post) + (1|study)",
            "Full: reduced + assembly covariates"),
  parameters = c(attr(logLik(reduced_fit), "df"), attr(logLik(fit), "df")),
  log_likelihood = c(as.numeric(logLik(reduced_fit)), as.numeric(logLik(fit))),
  AIC = c(AIC(reduced_fit), AIC(fit)),
  BIC = c(BIC(reduced_fit), BIC(fit)),
  convergence_code = c(reduced_fit$fit$convergence, fit$fit$convergence),
  positive_definite_hessian = c(isTRUE(reduced_fit$sdr$pdHess), isTRUE(fit$sdr$pdHess))
)
write.csv(fit_comparison,
          file.path(output_dir, "assembly_covariate_model_comparison.csv"), row.names = FALSE)

# Intermediate parsimonious model: retain non-host length only.
m_parsimonious <- glmmTMB(
  cbind(n_Em, n_Dis) ~ cohort + log(k_post) + log_non_host_length_z + (1 | study),
  ziformula = ~0,
  family = betabinomial(link = "logit"),
  data = pat,
  REML = FALSE
)
saveRDS(m_parsimonious, file.path(output_dir, "beta_binomial_model_parsimonious_non_host_length.rds"))
capture.output(
  summary(m_parsimonious),
  file = file.path(output_dir, "beta_binomial_model_parsimonious_summary_non_host_length.txt")
)
parsimonious_lrt <- anova(m_parsimonious, fit)
capture.output(
  parsimonious_lrt,
  file = file.path(output_dir, "parsimonious_vs_full_likelihood_ratio_test.txt")
)
parsimonious_comparison <- data.frame(
  model = c("Parsimonious: non_host_length only", "Full: all assembly covariates"),
  parameters = c(attr(logLik(m_parsimonious), "df"), attr(logLik(fit), "df")),
  log_likelihood = c(as.numeric(logLik(m_parsimonious)), as.numeric(logLik(fit))),
  AIC = c(AIC(m_parsimonious), AIC(fit)),
  BIC = c(BIC(m_parsimonious), BIC(fit)),
  convergence_code = c(m_parsimonious$fit$convergence, fit$fit$convergence),
  positive_definite_hessian = c(
    isTRUE(m_parsimonious$sdr$pdHess), isTRUE(fit$sdr$pdHess)
  )
)
write.csv(
  parsimonious_comparison,
  file.path(output_dir, "parsimonious_vs_full_model_comparison.csv"),
  row.names = FALSE
)

# Test whether the primary model generates fewer zero-emergence outcomes than observed.
set.seed(20260901)
simulated_residuals <- DHARMa::simulateResiduals(fittedModel = fit, n = 1000, plot = FALSE)
zero_test <- DHARMa::testZeroInflation(simulated_residuals, alternative = "greater")
capture.output(zero_test, file = file.path(output_dir, "dharma_zero_inflation_test.txt"))

zero_diagnostic <- data.frame(
  statistic = unname(zero_test$statistic),
  p_value = zero_test$p.value,
  alternative = zero_test$alternative,
  excess_zeros_detected = is.finite(zero_test$p.value) && zero_test$p.value < 0.05
)
write.csv(zero_diagnostic, file.path(output_dir, "dharma_zero_inflation_test.csv"), row.names = FALSE)

# Keep the non-zero-inflated fit primary. Fit ~1 only when DHARMa detects excess zeros.
model_comparison <- data.frame(
  model = "Primary: ziformula = ~0", parameters = attr(logLik(fit), "df"),
  log_likelihood = as.numeric(logLik(fit)), AIC = AIC(fit), BIC = BIC(fit)
)
if (zero_diagnostic$excess_zeros_detected) {
  fit_zi <- update(fit, ziformula = ~1)
  saveRDS(fit_zi, file.path(output_dir, "beta_binomial_model_zero_inflated_candidate.rds"))
  capture.output(
    summary(fit_zi),
    file = file.path(output_dir, "beta_binomial_model_zero_inflated_candidate_summary.txt")
  )
  model_comparison <- rbind(model_comparison, data.frame(
    model = "Candidate: ziformula = ~1", parameters = attr(logLik(fit_zi), "df"),
    log_likelihood = as.numeric(logLik(fit_zi)), AIC = AIC(fit_zi), BIC = BIC(fit_zi)
  ))
}
write.csv(model_comparison, file.path(output_dir, "zero_inflation_model_comparison.csv"), row.names = FALSE)

# Compute odds ratios for cohort comparisons using emmeans and contrasts, then save the results.
emm <- emmeans(fit, ~ cohort)
planned <- contrast(
  emm,
  method = list("MDRB vs rCDI" = c(-1, 1, 0), "Melanoma vs rCDI" = c(-1, 0, 1)),
  adjust = "holm"
)
or_table <- as.data.frame(summary(planned, infer = c(TRUE, TRUE), type = "response"))
write.csv(or_table, file.path(output_dir, "cohort_odds_ratios_holm.csv"), row.names = FALSE)

# Extract the study random-intercept standard deviation from the fitted model and save it.
variance_components <- VarCorr(fit)
study_stddev <- attr(variance_components$cond$study, "stddev")
if (is.null(study_stddev) || !"(Intercept)" %in% names(study_stddev)) {
  stop("Could not extract the study random-intercept SD from VarCorr(fit).")
}
study_sd <- data.frame(
  component = "study random-intercept SD",
  SD = unname(study_stddev["(Intercept)"])
)
write.csv(study_sd, file.path(output_dir, "study_random_intercept_sd.csv"), row.names = FALSE)

diagnostics <- data.frame(
  patients_exported = nrow(pat_all), patients_fitted = nrow(pat), studies = nlevels(pat$study),
  k_post_adjusted = TRUE,
  num_contigs_adjusted = TRUE,
  total_length_adjusted = TRUE,
  non_host_length_adjusted = TRUE,
    estimation = "maximum likelihood",
  zero_inflation = "none (ziformula = ~0)",
  log_likelihood = as.numeric(logLik(fit)),
  AIC = AIC(fit),
  BIC = BIC(fit),
  convergence_code = fit$fit$convergence,
  positive_definite_hessian = isTRUE(fit$sdr$pdHess)
)
write.csv(diagnostics, file.path(output_dir, "model_diagnostics.csv"), row.names = FALSE)

assembly_effect <- as.data.frame(summary(fit)$coefficients$cond)[
  c("log_num_contigs_z", "log_total_length_z", "log_non_host_length_z"), , drop = FALSE
]
assembly_effect$term <- c(
  "1-SD increase in log1p median recipient num_contigs",
  "1-SD increase in log1p median recipient total_length",
  "1-SD increase in log1p median recipient non_host_length"
)
assembly_effect$odds_ratio <- exp(assembly_effect$Estimate)
assembly_effect$OR_95CI_low <- exp(assembly_effect$Estimate - qnorm(0.975) * assembly_effect$`Std. Error`)
assembly_effect$OR_95CI_high <- exp(assembly_effect$Estimate + qnorm(0.975) * assembly_effect$`Std. Error`)
write.csv(assembly_effect, file.path(output_dir, "assembly_covariate_effects.csv"), row.names = FALSE)

print(or_table)
print(study_sd)
print(diagnostics)
print(assembly_effect)
