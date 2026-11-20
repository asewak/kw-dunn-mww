# main_plots.R -- 4 main presentation figures in base R (detailed ggplot2
# figures stay in outputs.R / appendix). Reads the long rejection-rate data frame
# `rej` (results/rejection_rates.rds). All rate figures use n3/n1 = 1 only.
#
#   main1_global_null.pdf   type I error vs n1, all groups N(0,1), incl. KW
#   main2_partial_null.pdf  type I error vs mu3 / sigma3, n1 = n2 = 50
#   main3_power.pdf         power (1 vs 2) vs mu3 / sigma3, delta = 0.5, n1 = n2 = 50
#   main4_dgps.pdf          densities of the location and scale group-3 families

# unibeCols may be missing (e.g. on the cluster): fall back to plain colours
# so master.R never fails just because of the palette.
if (requireNamespace("unibeCols", quietly = TRUE)) {
  library(unibeCols)
} else {
  message("unibeCols not installed: using fallback colours in main plots")
  unibeSaphire <- function() "#1F5FA6"
  unibeRed     <- function() "#E3001B"
  unibeApricot <- function() "#F29F05"
  unibeOcean   <- function() "#00A0B0"
  unibeGrey    <- function() "#7F7F7F"
  unibeGreenS  <- function() c("#CFE3D4", "#A3C9AE", "#6FAA80", "#3F8A58", "#1E6B3C")
}

# Same colour + shape per method in every figure.
main_methods <- c("Wilcoxon", "Dunn", "KW two-step", "KW + Wilcoxon", "KW")
main_labels  <- c(Wilcoxon = "MWW", Dunn = "Dunn", "KW two-step" = "KW-Dunn",
                  "KW + Wilcoxon" = "KW-MWW", KW = "KW")
main_cols    <- c(Wilcoxon = unibeSaphire(), Dunn = unibeRed(),
                  "KW two-step" = unibeApricot(),
                  "KW + Wilcoxon" = unibeOcean(), KW = unibeGrey())
main_pch     <- c(Wilcoxon = 16, Dunn = 17, "KW two-step" = 15,
                  "KW + Wilcoxon" = 18, KW = 4)

# One line + points + 95% MC error bars for one method.
draw_method <- function(d, m, xvar, jit = 1) {
  d <- d[d$method == m, ]
  d <- d[order(d[[xvar]]), ]
  if (!nrow(d)) return(invisible())
  x <- d[[xvar]] * jit
  lo <- pmax(0, d$rate - 1.96 * d$mcse); hi <- pmin(1, d$rate + 1.96 * d$mcse)
  segments(x, lo, x, hi, col = adjustcolor(main_cols[m], 0.6))
  lines(x, d$rate, col = main_cols[m], lwd = 1.8)
  points(x, d$rate, col = main_cols[m], pch = main_pch[m], cex = 1.1)
}

# Legend strip in the bottom outer margin. By default it is centred on the whole
# figure; with center_plot = TRUE it is centred on the plot region of the current
# (last drawn) panel instead.
draw_legend <- function(methods, cex = 1.1, center_plot = FALSE) {
  xc <- 0.5
  if (center_plot) {
    fig <- par("fig"); plt <- par("plt")
    xr  <- fig[1] + (fig[2] - fig[1]) * plt[1:2]   # plot region, device-normalised x
    xc  <- mean(xr)
  }
  op <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
            xaxs = "i", yaxs = "i", new = TRUE)
  on.exit(par(xaxs = op$xaxs, yaxs = op$yaxs), add = TRUE)
  plot.new()   # usr = 0..1 exactly, so coordinates are device-normalised
  # Drawn by hand (not legend()): legend() pads its box unevenly, so the visible
  # symbols + labels would sit off-centre. Here the visible content is centred on xc.
  labs <- main_labels[methods]
  cw   <- strwidth("0", cex = cex)          # one character width
  sym  <- 2.2 * cw                          # line + point sample
  gap  <- 0.5 * cw                          # sample <-> label
  sep  <- 1.6 * cw                          # between items
  tw   <- strwidth(labs, cex = cex)
  w    <- sym + gap + tw                    # width of each item
  x0   <- xc - (sum(w) + sep * (length(labs) - 1)) / 2
  y    <- 0.03
  for (i in seq_along(labs)) {
    m <- methods[i]
    segments(x0, y, x0 + sym, y, col = main_cols[m], lwd = 1.8)
    points(x0 + sym / 2, y, col = main_cols[m], pch = main_pch[m], cex = 1.1)
    text(x0 + sym + gap, y, labs[i], adj = c(0, 0.5), cex = cex)
    x0 <- x0 + w[i] + sep
  }
}

# Plot 1 ----------------------------------------------------------------------
main_plot1 <- function(rej, ratio = 1) {
  # location (mu3 = 0) and scale (sigma3 = 1) are the same global-null distribution
  g <- rej[rej$delta == 0 & rej$mu3 == 0 & rej$sigma3 == 1 & rej$ratio == ratio &
             rej$family != "skewed" & rej$method %in% main_methods, ]
  key <- interaction(g$n1, g$method, drop = TRUE)
  d <- do.call(rbind, lapply(split(g, key), function(x)
    data.frame(n1 = x$n1[1], method = x$method[1],
               rate = mean(x$rate), mcse = sqrt(mean(x$mcse^2) / nrow(x)))))
  ymax <- max(0.1, max(d$rate + 1.96 * d$mcse))
  par(mfrow = c(1, 1), mar = c(4.2, 4.6, 0.8, 1), oma = c(3, 0, 0, 0),
      mgp = c(2.6, 0.7, 0))
  ns <- sort(unique(d$n1))
  plot(NA, xlim = range(ns) * c(0.9, 1.1), ylim = c(0, ymax), log = "x",
       xaxt = "n", yaxt = "n",
       xlab = expression(n[1] == n[2] ~ "(log scale)"), ylab = "Type I error")
  axis(1, at = ns)
  axis(2, las = 2)
  abline(h = 0.05, lty = 2, col = "grey40")
  for (m in main_methods) draw_method(d, m, "n1")
  draw_legend(main_methods, cex = 0.85, center_plot = TRUE)
}

# Plots 2 and 3 ---------------------------------------------------------------
# one panel per group-3 family, n3 = n1; same layout for null and power.
main_rate_grid <- function(rej, delta, n1, ylab, title, null = TRUE, ratio = 1,
                           ref = NULL) {
  methods <- setdiff(main_methods, "KW")
  fams <- list(
    list(f = "location", x = "mu3",
         xlab = expression(mu[3] ~ "(group 3 location)"),
         lab = "Location"),
    list(f = "scale",    x = "sigma3",
         xlab = expression(sigma[3] ~ "(group 3 SD, log scale)"),
         lab = "Scale"))
  d <- rej[rej$delta == delta & rej$n1 == n1 & rej$ratio == ratio &
             rej$method %in% methods, ]
  # null: y axis in steps of 0.05, up to the next multiple covering the data
  ymax <- if (null) 0.05 * ceiling(max(0.1, max(d$rate + 1.96 * d$mcse)) / 0.05) else 1
  yat  <- if (null) seq(0, ymax, by = 0.05) else seq(0, 1, by = 0.2)
  par(mfrow = c(1, length(fams)), mar = c(4.2, 4.6, 0.8, 0.5), oma = c(3, 0, 0, 0),
      mgp = c(2.6, 0.7, 0))
  for (fm in fams) {
    dd <- d[d$family == fm$f, ]
    xs <- sort(unique(dd[[fm$x]]))
    plot(NA, xlim = range(xs), ylim = c(0, ymax), xaxt = "n", yaxt = "n",
         log = if (fm$x == "sigma3") "x" else "",
         xlab = fm$xlab, ylab = if (fm$f == "location") ylab else "")
    axis(1, at = xs, labels = xs)
    axis(2, at = yat, las = 2)
    if (null) abline(h = 0.05, lty = 2, col = "grey40")
    if (!is.null(ref)) abline(h = ref, lty = 2, col = "grey40")
    if (fm$x == "sigma3") abline(v = 1, lty = 3, col = "grey80")
    for (m in methods) draw_method(dd, m, fm$x)
  }
  draw_legend(methods)
}

main_plot2 <- function(rej, n1 = 50)
  main_rate_grid(rej, 0, n1, "Type I error",
                 sprintf("Partial null (group 1 = group 2), n1 = n2 = %d", n1))

main_plot3 <- function(rej, n1 = 50, delta = 0.5)
  main_rate_grid(rej, delta, n1, "Power (1 vs 2)",
                 sprintf("Power for 1 vs 2, group 2 ~ N(%s,1), n1 = n2 = %d", delta, n1),
                 null = FALSE,
                 ref = power.t.test(n = n1, delta = delta, sd = 1)$power)

# Plot 4: DGP densities --------------------------------------------------------
main_plot4 <- function(delta = 0.5) {
  # unibeGreenS shades, light -> dark (lightest shade dropped: too faint on white)
  g3col <- function(k) colorRampPalette(rev(unibeGreenS()[1:4]))(k)
  panel <- function(xlim, grp3, xlab, ylim = NULL) {
    # grp3: list(dens = function(x, par), pars, labs)
    ymax <- 0
    xs <- seq(xlim[1], xlim[2], length.out = 600)
    dens <- lapply(grp3$pars, function(p) grp3$dens(xs, p))
    ymax <- max(1.05 * dnorm(0), sapply(dens, max))
    if (is.null(ylim)) ylim <- c(0, 1.95 * ymax)
    plot(NA, xlim = xlim, ylim = ylim, xlab = xlab,
         ylab = "Density", yaxt = "n")
    axis(2, las = 2)
    cols <- g3col(length(dens))
    for (i in seq_along(dens)) lines(xs, dens[[i]], col = cols[i], lwd = 2)
    lines(xs, dnorm(xs), col = "black", lwd = 2.5)
    lines(xs, dnorm(xs, delta), col = "black", lwd = 2.5, lty = 2)
    legend("topright", bty = "n", cex = 0.7,
           legend = c(expression("Group 1: N(0,1)"),
                      as.expression(bquote("Group 2: N(" * .(delta) * ",1)")),
                      grp3$labs),
           col = c("black", "black", cols), lty = c(1, 2, rep(1, length(cols))),
           lwd = c(2.5, 2.5, rep(2, length(cols))))
  }
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 0.8, 0.5), oma = c(0, 0, 0, 0),
      mgp = c(2.6, 0.7, 0))
  mus <- c(0, 0.5, 1, 2, 3)
  panel(c(-4, 7),
        list(dens = function(x, p) dnorm(x, p), pars = mus,
             labs = as.expression(lapply(mus, function(v)
               bquote("Group 3: " * mu[3] == .(v))))), "")
  sg <- c(0.25, 0.5, 2, 4)
  panel(c(-10, 10),
        list(dens = function(x, p) dnorm(x, 0, p), pars = sg,
             labs = as.expression(lapply(sg, function(v)
               bquote("Group 3: " * sigma[3] == .(v))))), "", ylim = c(0, 2))
}

save_main_plots <- function(rej, cfg, n1 = 50, delta = 0.5) {
  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)
  f <- function(n) file.path(cfg$output_dir, n)
  pdf(f("main1_global_null.pdf"), width = 6, height = 4.8);   main_plot1(rej); dev.off()
  pdf(f("main2_partial_null.pdf"), width = 9, height = 4.5);     main_plot2(rej, n1); dev.off()
  pdf(f("main3_power.pdf"), width = 9, height = 4.5);            main_plot3(rej, n1, delta); dev.off()
  pdf(f("main4_dgps.pdf"), width = 9, height = 4.5);            main_plot4(delta); dev.off()
  invisible(NULL)
}
