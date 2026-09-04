#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
model_dir <- dirname(normalizePath(script_file))
input_dir <- file.path(model_dir, "output/model_2")
output_dir <- file.path(model_dir, "output/model_2")
transition_types <- c("persistence", "disappearance", "emergence", "transfer")

fits <- list()
results <- list()
study_sd <- list()
diagnostics <- list()
fit_statistics <- list()

for (transition_type in transition_types) {
  pat <- read.csv(file.path(input_dir, paste0(transition_type, ".csv")), stringsAsFactors = FALSE)
  required <- c("patient", "cohort", "study", "n_type", "n_total", "k_post")
  if (length(setdiff(required, names(pat)))) stop("Missing columns in ", transition_type, ".csv")
  if (anyNA(pat[, required]) || any(pat$n_type < 0) || any(pat$n_type > pat$n_total) ||
      any(pat$n_total <= 0) || any(pat$k_post <= 0)) stop("Invalid data for ", transition_type)

  pat$cohort <- factor(pat$cohort, levels = c("rCDI", "MDRB", "Melanoma"))
  pat$study <- factor(pat$study)
  fit <- glmmTMB(
    cbind(n_type, n_total - n_type) ~ cohort + log(k_post) + (1 | study),
    family = betabinomial(link = "logit"), data = pat
  )
  fits[[transition_type]] <- fit

  pearson_residuals <- residuals(fit, type = "pearson")
  pearson_chisq <- sum(pearson_residuals^2, na.rm = TRUE)
  residual_df <- df.residual(fit)
  fit_statistics[[transition_type]] <- data.frame(
    transition = transition_type,
    patients = nobs(fit),
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    BIC = BIC(fit),
    minus2_logLik = -2 * as.numeric(logLik(fit)),
    residual_df = residual_df,
    beta_binomial_dispersion = sigma(fit),
    pearson_chisq = pearson_chisq,
    pearson_chisq_df_ratio = pearson_chisq / residual_df,
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = isTRUE(fit$sdr$pdHess)
  )
  capture.output(summary(fit), file = file.path(output_dir, paste0(transition_type, "_model_summary.txt")))

  emm <- emmeans(fit, ~ cohort)
  contrasts <- contrast(emm, method = list(
    "MDRB vs rCDI" = c(-1, 1, 0),
    "Melanoma vs rCDI" = c(-1, 0, 1)
  ), adjust = "none")
  contrast_table <- as.data.frame(summary(contrasts, infer = c(TRUE, TRUE), type = "response"))
  contrast_table$transition <- transition_type
  results[[transition_type]] <- contrast_table

  vc <- VarCorr(fit)
  sd_value <- attr(vc$cond$study, "stddev")["(Intercept)"]
  study_sd[[transition_type]] <- data.frame(transition = transition_type, study_SD = unname(sd_value))
  diagnostics[[transition_type]] <- data.frame(
    transition = transition_type, patients = nrow(pat), studies = nlevels(pat$study),
    convergence_code = fit$fit$convergence,
    positive_definite_hessian = isTRUE(fit$sdr$pdHess)
  )
}

result_table <- do.call(rbind, results)
result_table$holm_p_across_four <- ave(result_table$p.value, result_table$contrast,
                                       FUN = function(x) p.adjust(x, method = "holm"))
result_table <- result_table[, c(
  "transition", "contrast", "odds.ratio", "asymp.LCL", "asymp.UCL",
  "p.value", "holm_p_across_four"
)]

write.csv(result_table, file.path(output_dir, "transition_model_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, study_sd), file.path(output_dir, "transition_study_sd.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(output_dir, "transition_model_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, fit_statistics), file.path(output_dir, "transition_model_fit_statistics.csv"), row.names = FALSE)
saveRDS(fits, file.path(output_dir, "transition_models.rds"))
print(result_table)
