# outputs.R -- extract rejection rates from the saved p-value objects and plot.

# Combine results (list over n1) into one long data frame of rejection rates.
extract_rejections <- function(results, cfg) {
  a <- cfg$alpha
  rows <- lapply(results, function(res) {
    d <- res$design
    do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
      p <- res$pvals[[i]]
      rate <- c(KW = mean(p[, "kw"] < a),
                "KW two-step" = mean(p[, "kw"] < a & p[, "dunn"] < a),
                Wilcoxon = mean(p[, "wilcox"] < a),
                Dunn = mean(p[, "dunn"] < a))
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
  library(ggplot2)
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
    ggplot2::ggsave(f, plot_rejections(rej, fam, dl), width = 12, height = 9)
  }
  invisible(NULL)
}
