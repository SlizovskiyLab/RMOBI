#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file))
output_dir <- file.path(script_dir, "output", "leave_one_study_out")
pat <- read.csv(file.path(output_dir, "loso_glmm_patients.csv"), stringsAsFactors = FALSE)
pat$cohort <- factor(pat$cohort, levels = c("rCDI", "MDRB", "Melanoma"))
pat$study <- factor(pat$study)

fit_one <- function(data, excluded_study) {
  data$study <- droplevels(data$study)
  fit <- glmmTMB(
    cbind(n_Em, n_Dis) ~ cohort + log(k_post) + (1 | study),
    family = betabinomial(link = "logit"), data = data
  )
  emm <- emmeans(fit, ~ cohort)
  contrasts <- contrast(emm, method = list(
    "MDRB vs rCDI" = c(-1, 1, 0),
    "Melanoma vs rCDI" = c(-1, 0, 1)
  ), adjust = "none")
  result <- as.data.frame(summary(contrasts, infer = c(TRUE, TRUE), type = "response"))
  link_result <- as.data.frame(summary(contrasts, infer = c(TRUE, TRUE), type = "link"))
  result$log_odds_ratio <- link_result$estimate
  result$excluded_study <- excluded_study
  result$patients <- nrow(data)
  result$studies <- nlevels(data$study)
  result$study_SD <- unname(attr(VarCorr(fit)$cond$study, "stddev")["(Intercept)"])
  result$dispersion <- sigma(fit)
  result$convergence_code <- fit$fit$convergence
  result$positive_definite_hessian <- isTRUE(fit$sdr$pdHess)
  list(fit = fit, result = result)
}

fits <- list()
rows <- list()
full <- fit_one(pat, "None (full model)")
fits[["full_model"]] <- full$fit
rows[["full_model"]] <- full$result
for (excluded in sort(unique(as.character(pat$study)))) {
  current <- fit_one(subset(pat, study != excluded), excluded)
  fits[[excluded]] <- current$fit
  rows[[excluded]] <- current$result
}

results <- do.call(rbind, rows)
results <- results[, c(
  "excluded_study", "contrast", "log_odds_ratio", "odds.ratio", "SE",
  "asymp.LCL", "asymp.UCL", "z.ratio", "p.value", "patients", "studies",
  "study_SD", "dispersion", "convergence_code", "positive_definite_hessian"
)]
write.csv(results, file.path(output_dir, "loso_glmm_coefficients.csv"), row.names = FALSE)
saveRDS(fits, file.path(output_dir, "loso_glmm_models.rds"))

labels <- c("None (full model)", sort(unique(as.character(pat$study))))
keys <- as.vector(rbind(paste(labels, "MDRB vs rCDI"), paste(labels, "Melanoma vs rCDI")))
plot_data <- results[match(keys, paste(results$excluded_study, results$contrast)), ]
colors <- c("MDRB vs rCDI" = "#9f2b23", "Melanoma vs rCDI" = "#555555")
y_base <- rev(seq_along(labels))
y <- rep(y_base, each = 2) + rep(c(0.13, -0.13), length(labels))
x_limits <- range(c(plot_data$asymp.LCL, plot_data$asymp.UCL, 1), finite = TRUE)

draw_forest <- function() {
  par(mar = c(4.5, 10, 3.5, 1))
  plot(plot_data$odds.ratio, y, log = "x", xlim = x_limits,
       ylim = c(0.5, length(labels) + 0.5), yaxt = "n", pch = 19,
       col = colors[plot_data$contrast], xlab = "Cohort odds ratio (95% CI)", ylab = "")
  segments(plot_data$asymp.LCL, y, plot_data$asymp.UCL, y,
           col = colors[plot_data$contrast], lwd = 1.5)
  axis(2, at = y_base, labels = labels, las = 1, cex.axis = 0.72)
  abline(v = 1, lty = 2, col = "grey45")
  legend("top", inset = c(0, -0.12), legend = names(colors), col = colors,
         pch = 19, horiz = TRUE, xpd = NA, bty = "n")
}

png(file.path(output_dir, "loso_glmm_coefficient_forest.png"),
    width = 2100, height = 1800, res = 300)
draw_forest(); dev.off()
svg(file.path(output_dir, "loso_glmm_coefficient_forest.svg"), width = 7, height = 6)
draw_forest(); dev.off()
print(results)
