# outputs.R -- extract rejection rates from the saved p-value objects and plot.

library(ggplot2)

# Combine results (list over n1) into one long data frame of rejection rates.
extract_rejections <- function(results, cfg) {
  a <- cfg$alpha
  rows <- lapply(results, function(res) {
    d <- res$design
    do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
      p <- res$pvals[[i]]
      rate <- c(KW = mean(p[, "kw"] < a),
                "KW two-step" = mean(p[, "kw"] < a & p[, "dunn"] < a),
                "KW + Wilcoxon" = mean(p[, "kw"] < a & p[, "wilcox"] < a),
                Wilcoxon = mean(p[, "wilcox"] < a),
                Dunn = mean(p[, "dunn"] < a),
                "Pairwise t" = mean(p[, "pairt"] < a),
                "Pairwise t (no pool)" = mean(p[, "pairt_np"] < a))
      # KW alone is only meaningful under the global null
      global_null <- d$delta[i] == 0 && d$family[i] != "skewed" &&
        d$mu3[i] == 0 && d$sigma3[i] == 1
      if (!global_null) rate["KW"] <- NA
      data.frame(d[i, ], method = names(rate), rate = unname(rate),
                 mcse = sqrt(unname(rate) * (1 - unname(rate)) / nrow(p)),
                 row.names = NULL)
    }))
  })
  out <- do.call(rbind, rows)
  out[!is.na(out$rate), ]
}

# One figure per family x delta: rejection rate vs mu3 (location/skewed) or
# sigma3 (scale); columns = n1, rows = n3/n1.
plot_rejections <- function(rej, family, delta) {
  d <- rej[rej$family == family & rej$delta == delta, ]
  xvar <- if (family == "scale") "sigma3" else "mu3"
  d$x <- d[[xvar]]
  d$n1f <- factor(d$n1, labels = paste0("n1 = n2 = ", sort(unique(d$n1))))
  d$ratiof <- factor(d$ratio, labels = paste0("n3/n1 = ", sort(unique(d$ratio))))
  ttl <- sprintf("%s group 3, shift of group 2 = %s (%s)", family, delta,
                 if (delta == 0) "null for 1 vs 2" else "power for 1 vs 2")
  p <- ggplot(d, aes(x, rate, colour = method, shape = method)) +
    geom_linerange(aes(ymin = pmax(0, rate - 1.96 * mcse),
                       ymax = pmin(1, rate + 1.96 * mcse)), alpha = 0.4) +
    geom_line() + geom_point() +
    facet_grid(ratiof ~ n1f) +
    labs(x = xvar, y = "Rejection rate", title = ttl) +
    theme_bw() + theme(legend.position = "bottom")
  if (family == "scale") p <- p + scale_x_continuous(trans = "log2")
  if (delta == 0) {
    # nominal level reference line (only meaningful for null settings)
    p + geom_hline(yintercept = 0.05, linetype = 2, colour = "grey50")
  } else {
    p + coord_cartesian(ylim = c(0, 1))
  }
}

# Rejection rates (with Monte Carlo SEs) as a plain CSV.
save_rejection_csv <- function(rej, cfg) {
  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(rej, file.path(cfg$output_dir, "rejection_rates.csv"), row.names = FALSE)
}

# Record how the results were produced: settings, R/package versions, date,
# git commit.
save_provenance <- function(cfg) {
  commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE,
                            ignore.stderr = TRUE), error = function(e) NA,
                     warning = function(w) NA)
  prov <- list(cfg = cfg, date = Sys.time(), git_commit = commit,
               dunn.test = as.character(packageVersion("dunn.test")),
               session = sessionInfo())
  dir.create(cfg$results_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(prov, file.path(cfg$results_dir, "provenance.rds"))
  writeLines(c(paste("date:", format(prov$date)), paste("git:", commit),
               paste("dunn.test:", prov$dunn.test),
               capture.output(print(prov$session))),
             file.path(cfg$results_dir, "provenance.txt"))
}

save_all_plots <- function(rej, cfg) {
  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)
  for (fam in unique(rej$family)) for (dl in cfg$delta_grid) {
    f <- file.path(cfg$output_dir, sprintf("rejection_%s_delta%s.pdf", fam, dl))
    ggsave(f, plot_rejections(rej, fam, dl), width = 12, height = 9)
  }
  invisible(NULL)
}

# ---- summary figures (7 plots) -----------------------------------------------
# Consistent method colours / order across all summary plots.
method_cols <- c(Wilcoxon = "#0072B2", Dunn = "#D55E00",
                 "KW two-step" = "#009E73", "KW + Wilcoxon" = "#CC79A7",
                 "Pairwise t" = "#E69F00", "Pairwise t (no pool)" = "#56B4E9",
                 KW = "#7F7F7F")

summary_ratios <- c(0.5, 1, 2)   # n3/n1 shown in the summary plots

summary_prep <- function(d) {
  d <- d[d$ratio %in% summary_ratios, ]
  d$method <- factor(d$method, levels = names(method_cols))
  d$n1f    <- factor(d$n1, levels = sort(unique(d$n1)),
                     labels = paste0("n1 = n2 = ", sort(unique(d$n1))))
  d$ratiof <- factor(d$ratio, levels = sort(unique(d$ratio)),
                     labels = paste0("n3/n1 = ", sort(unique(d$ratio))))
  d
}

summary_base <- function(d, xvar, xlab, ylab, ref = TRUE) {
  p <- ggplot(d, aes(.data[[xvar]], rate, colour = method, shape = method)) +
    geom_linerange(aes(ymin = pmax(0, rate - 1.96 * mcse),
                       ymax = pmin(1, rate + 1.96 * mcse)), alpha = 0.7) +
    geom_line() + geom_point(size = 1.8) +
    scale_colour_manual(values = method_cols, drop = TRUE) +
    labs(x = xlab, y = ylab, colour = NULL, shape = NULL) +
    theme_bw() + theme(legend.position = "bottom")
  if (ref) p <- p + geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40")
  p
}

# Plot 1: global null; location (mu3 = 0) and scale (sigma3 = 1) runs are the
# same distribution, so pool them.
plot_summary1 <- function(rej) {
  g <- rej[rej$delta == 0 & rej$mu3 == 0 & rej$sigma3 == 1 &
             rej$family != "skewed", ]
  key <- interaction(g$n1, g$ratio, g$method)
  pooled <- do.call(rbind, lapply(split(g, key), function(x)
    data.frame(x[1, c("n1", "ratio", "method")], rate = mean(x$rate),
               mcse = sqrt(mean(x$mcse^2) / nrow(x)))))
  d <- summary_prep(pooled)
  summary_base(d, "n1", "n1 = n2 (log scale)", "Type I error") +
    scale_x_continuous(trans = "log10", breaks = sort(unique(d$n1))) +
    facet_grid(. ~ ratiof) +
    ggtitle("Global null: all groups N(0,1)")
}

# Plots 2 and 3 (and 4, below): partial nulls (delta = 0)
plot_summary2 <- function(rej) {
  d <- summary_prep(rej[rej$delta == 0 & rej$family == "location" &
                          rej$method != "KW", ])
  summary_base(d, "mu3", "mu3 (location of group 3)", "Type I error") +
    facet_grid(ratiof ~ n1f) +
    ggtitle("Partial null, location: group 1 = group 2, group 3 ~ N(mu3, 1)")
}

plot_summary3 <- function(rej) {
  d <- summary_prep(rej[rej$delta == 0 & rej$family == "scale" &
                          rej$method != "KW", ])
  summary_base(d, "sigma3", "sigma3 (SD of group 3, log scale)", "Type I error") +
    scale_x_continuous(trans = "log2", breaks = sort(unique(d$sigma3))) +
    facet_grid(ratiof ~ n1f) +
    ggtitle("Partial null, scale: group 1 = group 2, group 3 ~ N(0, sigma3^2)")
}

# Plot 5: power vs mu3, location family, one n1
plot_summary5 <- function(rej, n1 = 50) {
  d <- rej[rej$delta > 0 & rej$family == "location" & rej$n1 == n1 &
             rej$method != "KW", ]
  d <- summary_prep(d)
  d$deltaf <- factor(d$delta, labels = paste0("delta = ", sort(unique(d$delta))))
  summary_base(d, "mu3", "mu3 (location of group 3)", "Power (1 vs 2)", ref = FALSE) +
    facet_grid(ratiof ~ deltaf) + coord_cartesian(ylim = c(0, 1)) +
    ggtitle(sprintf("Power, location family, n1 = n2 = %d", n1))
}

# Plot 6: power vs sigma3, scale family, one n1
plot_summary6 <- function(rej, n1 = 50) {
  d <- rej[rej$delta > 0 & rej$family == "scale" & rej$n1 == n1 &
             rej$method != "KW", ]
  d <- summary_prep(d)
  d$deltaf <- factor(d$delta, labels = paste0("delta = ", sort(unique(d$delta))))
  summary_base(d, "sigma3", "sigma3 (SD of group 3, log scale)", "Power (1 vs 2)",
               ref = FALSE) +
    scale_x_continuous(trans = "log2", breaks = sort(unique(d$sigma3))) +
    facet_grid(ratiof ~ deltaf) + coord_cartesian(ylim = c(0, 1)) +
    ggtitle(sprintf("Power, scale family, n1 = n2 = %d", n1))
}

# Plot 4: partial null, skewed group 3 (mu3 + Exp(1) - 1); delta = 0
plot_summary4 <- function(rej) {
  d <- summary_prep(rej[rej$delta == 0 & rej$family == "skewed" &
                          rej$method != "KW", ])
  summary_base(d, "mu3", "mu3 (location of group 3)", "Type I error") +
    facet_grid(ratiof ~ n1f) +
    ggtitle("Partial null, skewed: group 1 = group 2, group 3 ~ mu3 + (Exp(1) - 1)")
}

# Plot 7: power vs mu3, skewed family, one n1
plot_summary7 <- function(rej, n1 = 50) {
  d <- rej[rej$delta > 0 & rej$family == "skewed" & rej$n1 == n1 &
             rej$method != "KW", ]
  d <- summary_prep(d)
  d$deltaf <- factor(d$delta, labels = paste0("delta = ", sort(unique(d$delta))))
  summary_base(d, "mu3", "mu3 (location of group 3)", "Power (1 vs 2)", ref = FALSE) +
    facet_grid(ratiof ~ deltaf) + coord_cartesian(ylim = c(0, 1)) +
    ggtitle(sprintf("Power, skewed family, n1 = n2 = %d", n1))
}

save_summary_plots <- function(rej, cfg) {
  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)
  f <- function(n) file.path(cfg$output_dir, n)
  ggsave(f("sum1_global_null.pdf"),     plot_summary1(rej), width = 12, height = 4)
  ggsave(f("sum2_partial_null_location.pdf"), plot_summary2(rej), width = 14, height = 8)
  ggsave(f("sum3_partial_null_scale.pdf"),    plot_summary3(rej), width = 14, height = 8)
  ggsave(f("sum4_partial_null_skewed.pdf"),   plot_summary4(rej), width = 14, height = 8)
  ggsave(f("sum5_power_location.pdf"),  plot_summary5(rej), width = 10, height = 8)
  ggsave(f("sum6_power_scale.pdf"),     plot_summary6(rej), width = 10, height = 8)
  ggsave(f("sum7_power_skewed.pdf"),          plot_summary7(rej), width = 10, height = 8)
  invisible(NULL)
}
