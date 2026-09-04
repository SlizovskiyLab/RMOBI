#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file))
output_dir <- file.path(script_dir, "output", "model_3")
pat <- read.csv(file.path(output_dir, "patient_mge_emergence_disappearance.csv"), stringsAsFactors = FALSE)

assembly_columns <- c(
  "median_recipient_num_contigs", "median_recipient_total_length",
  "median_recipient_non_host_length"
)
required <- c("patient", "study", "cohort", "MGE_class", "n_Em", "n_Dis", "k_post",
              assembly_columns)
if (length(setdiff(required, names(pat))) || anyNA(pat[, required])) stop("Invalid Model 4 input.")
if (any(pat$n_Em < 0 | pat$n_Dis < 0 | pat$n_Em + pat$n_Dis <= 0 | pat$k_post <= 0)) {
  stop("Model 4 counts or k_post are invalid.")
}

# ICE is excluded from all Model 4 fits and reported separately for audit.
excluded_ice <- pat[toupper(trimws(pat$MGE_class)) == "ICE", , drop = FALSE]
pat <- pat[toupper(trimws(pat$MGE_class)) != "ICE", , drop = FALSE]
write.csv(excluded_ice, file.path(output_dir, "excluded_ice_rows.csv"), row.names = FALSE)
if (!nrow(pat)) stop("No Model 4 rows remain after excluding ICE.")

pat$cohort <- factor(pat$cohort, levels = c("rCDI", "MDRB", "Melanoma"))
pat$MGE_class <- relevel(factor(pat$MGE_class), ref = "plasmid")
pat$study <- factor(pat$study)
pat$patient <- factor(pat$patient)

# Standardize assembly metrics once per patient so patients represented in more
# MGE classes do not receive extra weight when calculating means and SDs.
patient_assembly <- unique(pat[c("patient", assembly_columns)])
for (column in assembly_columns) {
  patient_assembly[[column]] <- as.numeric(patient_assembly[[column]])
  if (any(!is.finite(patient_assembly[[column]])) || any(patient_assembly[[column]] <= 0)) {
    stop(column, " must contain positive finite values.")
  }
}
patient_assembly$log_num_contigs_z <- as.numeric(
  scale(log1p(patient_assembly$median_recipient_num_contigs))
)
patient_assembly$log_total_length_z <- as.numeric(
  scale(log1p(patient_assembly$median_recipient_total_length))
)
patient_assembly$log_non_host_length_z <- as.numeric(
  scale(log1p(patient_assembly$median_recipient_non_host_length))
)
z_columns <- c("log_num_contigs_z", "log_total_length_z", "log_non_host_length_z")
for (column in z_columns) {
  pat[[column]] <- patient_assembly[[column]][match(pat$patient, patient_assembly$patient)]
}

reduced <- glmmTMB(
  cbind(n_Em, n_Dis) ~ cohort * MGE_class + log(k_post) +
    (1 | study) + (1 | patient),
  ziformula = ~0, family = betabinomial(link = "logit"), data = pat, REML = FALSE
)
parsimonious <- glmmTMB(
  cbind(n_Em, n_Dis) ~ cohort * MGE_class + log(k_post) +
    log_non_host_length_z + (1 | study) + (1 | patient),
  ziformula = ~0, family = betabinomial(link = "logit"), data = pat, REML = FALSE
)
full <- glmmTMB(
  cbind(n_Em, n_Dis) ~ cohort * MGE_class + log(k_post) +
    log_num_contigs_z + log_total_length_z + log_non_host_length_z +
    (1 | study) + (1 | patient),
  ziformula = ~0, family = betabinomial(link = "logit"), data = pat, REML = FALSE
)

# The primary interaction test compares models with identical technical covariates.
interaction_null <- glmmTMB(
  cbind(n_Em, n_Dis) ~ cohort + MGE_class + log(k_post) +
    log_num_contigs_z + log_total_length_z + log_non_host_length_z +
    (1 | study) + (1 | patient),
  ziformula = ~0, family = betabinomial(link = "logit"), data = pat, REML = FALSE
)
interaction_test <- as.data.frame(anova(interaction_null, full))
interaction_test$model <- rownames(interaction_test)
rownames(interaction_test) <- NULL
interaction_test$hypothesis <- c(NA, "cohort × MGE class interaction = 0")

candidate_comparison <- as.data.frame(anova(reduced, parsimonious, full))
candidate_comparison$model <- rownames(candidate_comparison)
rownames(candidate_comparison) <- NULL
candidate_comparison$comparison <- c(
  "Reference", "Parsimonious vs reduced", "Full vs parsimonious"
)
cell_estimates <- as.data.frame(summary(
  emmeans(full, ~ cohort * MGE_class), infer = c(TRUE, FALSE), type = "response"
))
posthoc_emm <- emmeans(full, ~ cohort | MGE_class)
posthoc <- contrast(posthoc_emm, method = list(
  "MDRB vs rCDI" = c(-1, 1, 0),
  "Melanoma vs rCDI" = c(-1, 0, 1)
), adjust = "none")
posthoc_results <- as.data.frame(summary(
  posthoc, infer = c(TRUE, TRUE), type = "response", adjust = "none"
))
posthoc_results$holm_p_across_8 <- p.adjust(posthoc_results$p.value, method = "holm")
posthoc_results$significant_holm_0.05 <- posthoc_results$holm_p_across_8 < 0.05
fixed_effects <- as.data.frame(summary(full)$coefficients$cond)
fixed_effects$term <- rownames(fixed_effects)
rownames(fixed_effects) <- NULL

vc <- VarCorr(full)$cond
random_effects <- data.frame(
  group = c("study", "patient"),
  SD = c(unname(attr(vc$study, "stddev")["(Intercept)"]),
         unname(attr(vc$patient, "stddev")["(Intercept)"]))
)
diagnostics <- data.frame(
  rows = nrow(pat), patients = nlevels(pat$patient), studies = nlevels(pat$study),
  MGE_classes = nlevels(pat$MGE_class), logLik = as.numeric(logLik(full)),
  AIC = AIC(full), BIC = BIC(full), dispersion = sigma(full),
  convergence_code = full$fit$convergence,
  positive_definite_hessian = isTRUE(full$sdr$pdHess)
)

capture.output({
  cat("PRIMARY HYPOTHESIS TEST\n"); print(interaction_test, row.names = FALSE)
  cat("\nFULL MODEL\n"); print(summary(full))
  cat("\nPARSIMONIOUS MODEL\n"); print(summary(parsimonious))
  cat("\nREDUCED MODEL\n"); print(summary(reduced))
  cat("\nTHREE-MODEL ASSEMBLY-COVARIATE COMPARISON\n")
  print(candidate_comparison, row.names = FALSE)
  cat("\nRANDOM-EFFECT SDs\n"); print(random_effects, row.names = FALSE)
  cat("\nPOST-HOC COHORT COMPARISONS WITHIN MGE CLASS\n")
  print(posthoc_results, row.names = FALSE)
  cat("\nDIAGNOSTICS\n"); print(diagnostics, row.names = FALSE)
}, file = file.path(output_dir, "mge_class_interaction_model_summary.txt"))

write.csv(interaction_test, file.path(output_dir, "interaction_likelihood_ratio_test.csv"), row.names = FALSE)
write.csv(candidate_comparison,
          file.path(output_dir, "assembly_covariate_model_comparison.csv"), row.names = FALSE)
write.csv(fixed_effects, file.path(output_dir, "full_model_fixed_effects.csv"), row.names = FALSE)
write.csv(cell_estimates, file.path(output_dir, "cohort_mge_class_estimated_probabilities.csv"), row.names = FALSE)
write.csv(posthoc_results, file.path(output_dir, "within_mge_class_cohort_contrasts_holm.csv"), row.names = FALSE)
write.csv(random_effects, file.path(output_dir, "random_effect_sds.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_dir, "model_diagnostics.csv"), row.names = FALSE)
saveRDS(
  list(full = full, parsimonious = parsimonious, reduced = reduced,
       interaction_null = interaction_null),
  file.path(output_dir, "mge_class_interaction_models.rds")
)

plot_classes <- c("plasmid", "likely IS/TE", "prophage", "virus")
plot_cohorts <- c("rCDI", "MDRB", "Melanoma")
plot_data <- subset(cell_estimates, MGE_class %in% plot_classes)
plot_data$x <- match(plot_data$MGE_class, plot_classes) +
  c(rCDI = -0.18, MDRB = 0, Melanoma = 0.18)[as.character(plot_data$cohort)]
cohort_colors <- c(rCDI = "#174f88", MDRB = "#9f2b23", Melanoma = "#777777")

draw_probability_plot <- function() {
  par(mar = c(5, 5, 3.5, 1))
  plot(NA, xlim = c(0.55, 4.45), ylim = c(0, 1), xaxt = "n",
       xlab = "MGE class", ylab = "Adjusted probability of emergence")
  axis(1, at = seq_along(plot_classes),
       labels = c("Plasmid", "Likely IS/TE", "Prophage", "Virus"))
  segments(plot_data$x, plot_data$asymp.LCL, plot_data$x, plot_data$asymp.UCL,
           col = cohort_colors[plot_data$cohort], lwd = 1.6)
  segments(plot_data$x - 0.035, plot_data$asymp.LCL,
           plot_data$x + 0.035, plot_data$asymp.LCL,
           col = cohort_colors[plot_data$cohort], lwd = 1.6)
  segments(plot_data$x - 0.035, plot_data$asymp.UCL,
           plot_data$x + 0.035, plot_data$asymp.UCL,
           col = cohort_colors[plot_data$cohort], lwd = 1.6)
  points(plot_data$x, plot_data$prob, pch = 19,
         col = cohort_colors[plot_data$cohort], cex = 1.15)
  legend("top", inset = c(0, -0.16), legend = plot_cohorts,
         col = cohort_colors[plot_cohorts], pch = 19, horiz = TRUE, xpd = NA, bty = "n")
}

png(file.path(output_dir, "adjusted_emergence_probability_by_cohort_mge.png"),
    width = 1800, height = 1300, res = 300)
draw_probability_plot(); dev.off()
svg(file.path(output_dir, "adjusted_emergence_probability_by_cohort_mge.svg"),
    width = 6, height = 4.33)
draw_probability_plot(); dev.off()
print(interaction_test)
print(diagnostics)
