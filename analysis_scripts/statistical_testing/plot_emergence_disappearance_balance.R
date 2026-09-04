#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(emmeans))

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file))
output_dir <- file.path(script_dir, "output", "model_1")
fit <- readRDS(file.path(output_dir, "beta_binomial_model.rds"))
pat <- read.csv(file.path(output_dir, "emergence_disappearance_patients.csv"))
pat$log_num_contigs_z <- as.numeric(scale(log1p(pat$median_recipient_num_contigs)))
pat$log_total_length_z <- as.numeric(scale(log1p(pat$median_recipient_total_length)))
pat$log_non_host_length_z <- as.numeric(scale(log1p(pat$median_recipient_non_host_length)))
pat <- subset(pat, n_Em + n_Dis > 0)
pat$cohort <- factor(pat$cohort, levels = c("rCDI", "MDRB", "Melanoma"))
pat$study <- factor(pat$study)

# Cohort-specific adjusted emergence/disappearance odds from Model 1.
cohort_emm <- as.data.frame(summary(emmeans(fit, ~ cohort, data = pat), infer = c(TRUE, FALSE)))
cohort_emm$odds <- exp(cohort_emm$emmean)
cohort_emm$lower <- exp(cohort_emm$asymp.LCL)
cohort_emm$upper <- exp(cohort_emm$asymp.UCL)
cohort_emm <- cohort_emm[match(c("rCDI", "Melanoma", "MDRB"), cohort_emm$cohort), ]

## PNG (high-res raster)
png(file.path(output_dir, "cohort_emergence_disappearance_odds_forest.png"),
    width = 1800, height = 1200, res = 300)
par(mar = c(4.5, 7, 1, 1))
y <- seq_len(nrow(cohort_emm))
plot(cohort_emm$odds, y, log = "x", xlim = range(c(cohort_emm$lower, cohort_emm$upper, 1)),
     ylim = c(0.5, nrow(cohort_emm) + 0.5), pch = 19, yaxt = "n",
     xlab = "Emergence-to-disappearance odds (95% CI)", ylab = "")
segments(cohort_emm$lower, y, cohort_emm$upper, y, lwd = 2)
axis(2, at = y, labels = cohort_emm$cohort, las = 1)
abline(v = 1, lty = 2, col = "grey40")
dev.off()

## SVG (vector) - use same physical size as PNG (pixels/res -> inches)
svg(file.path(output_dir, "cohort_emergence_disappearance_odds_forest.svg"),
    width = 1800/300, height = 1200/300)
par(mar = c(4.5, 7, 1, 1))
y <- seq_len(nrow(cohort_emm))
plot(cohort_emm$odds, y, log = "x", xlim = range(c(cohort_emm$lower, cohort_emm$upper, 1)),
     ylim = c(0.5, nrow(cohort_emm) + 0.5), pch = 19, yaxt = "n",
     xlab = "Emergence-to-disappearance odds (95% CI)", ylab = "")
segments(cohort_emm$lower, y, cohort_emm$upper, y, lwd = 2)
axis(2, at = y, labels = cohort_emm$cohort, las = 1)
abline(v = 1, lty = 2, col = "grey40")
dev.off()

# Patient-level balance; the fitted relationship uses the model's log(k_post) scale.
pat$delta_p <- log((pat$n_Em + 0.5) / (pat$n_Dis + 0.5))
colors <- c(rCDI = "#174f88", MDRB = "#8f241e", Melanoma = "#8a8787")
line_fit <- lm(delta_p ~ log(k_post), data = pat)
x_line <- seq(min(pat$k_post), max(pat$k_post), length.out = 200)

## PNG (high-res raster)
png(file.path(output_dir, "delta_vs_post_fmt_samples.png"),
    width = 1800, height = 1200, res = 300)
set.seed(20260825)
par(mar = c(4.5, 6, 3.5, 1))
plot(jitter(pat$k_post, amount = 0.06), pat$delta_p, pch = 19,
     col = colors[pat$cohort], xaxt = "n", xlim = c(0.8, 4.2),
     xlab = "Number of post-FMT samples",
     ylab = "delta_p = log[(n_Em + 0.5)/(n_Dis + 0.5)]")
axis(1, at = sort(unique(pat$k_post)))
abline(h = 0, lty = 2, col = "grey60")
lines(x_line, predict(line_fit, newdata = data.frame(k_post = x_line)), lwd = 2)
legend("top", inset = c(0, -0.18), legend = names(colors), col = colors,
       pch = 19, horiz = TRUE, xpd = NA, bty = "n")
dev.off()

## SVG (vector) - use same physical size as PNG (pixels/res -> inches)
svg(file.path(output_dir, "delta_vs_post_fmt_samples.svg"),
    width = 1800/300, height = 1200/300)
set.seed(20260825)
par(mar = c(4.5, 6, 3.5, 1))
plot(jitter(pat$k_post, amount = 0.06), pat$delta_p, pch = 19,
     col = colors[pat$cohort], xaxt = "n", xlim = c(0.8, 4.2),
     xlab = "Number of post-FMT samples",
     ylab = "delta_p = log[(n_Em + 0.5)/(n_Dis + 0.5)]")
axis(1, at = sort(unique(pat$k_post)))
abline(h = 0, lty = 2, col = "grey60")
lines(x_line, predict(line_fit, newdata = data.frame(k_post = x_line)), lwd = 2)
legend("top", inset = c(0, -0.18), legend = names(colors), col = colors,
       pch = 19, horiz = TRUE, xpd = NA, bty = "n")
dev.off()
