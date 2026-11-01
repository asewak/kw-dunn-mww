# methods.R -- the tests. Each replicate is run through all methods and the
# p-values are saved in an object; rejections are extracted later (outputs.R).
#
# Methods:
#   kw        stats::kruskal.test on all three groups (reported only in global null)
#   kw2step   KW rejects AND Dunn 1 vs 2 rejects (formed in outputs.R)
#   wilcox    stats::wilcox.test, 1 vs 2 (ranks from the pair only)
#   dunn      dunn.test::dunn.test, 1 vs 2 (ranks pooled over all three groups)

# Dunn's test from the dunn.test package. We take the Z statistics, form
# two-sided p-values ourselves (dunn.test returns one-sided P by default) and
# adjust across the 3 pairs. Element 1 of the result is the 1 vs 2 comparison.
dunn_pvals <- function(y, g, p_adjust = "none") {
  # dunn.test prints a banner (stdout and/or messages); silence both
  invisible(utils::capture.output(suppressMessages(
    d <- dunn.test::dunn.test(y, g, method = "none", kw = FALSE, table = FALSE))))
  # put pairs in a fixed order (1-2, 1-3, 2-3) whatever labels/order the
  # package uses, so element 1 is always the 1 vs 2 comparison
  lab <- vapply(strsplit(d$comparisons, " - ", fixed = TRUE),
                function(s) paste(sort(trimws(s)), collapse = "-"), "")
  ord <- match(c("1-2", "1-3", "2-3"), lab)
  stopifnot(!anyNA(ord))
  p.adjust(2 * pnorm(-abs(d$Z[ord])), method = p_adjust)
}

# All methods on one replicate; returns named vector of p-values.
run_methods <- function(y, n, p_adjust = "none") {
  g <- rep(1:3, times = n)
  c(kw     = kruskal.test(y, g)$p.value,
    wilcox = wilcox.test(y[g == 1], y[g == 2])$p.value,
    dunn   = dunn_pvals(y, g, p_adjust)[1])
}

# Run all replicates of one scenario. Returns a reps x 3 matrix of p-values.
run_scenario <- function(datalist, n, cfg) {
  # forking costs ~0.2s per call, which swamps the work at small reps, so
  # only use multiple cores when there are enough replications to pay for it.
  # This is the only mclapply in the pipeline (never nested).
  cores <- if (length(datalist) >= cfg$parallel_min_reps) cfg$mc_cores else 1
  res <- parallel::mclapply(datalist, run_methods, n = n,
                            p_adjust = cfg$dunn_p_adjust,
                            mc.cores = cores)
  # mclapply returns try-error objects (or NULL) instead of stopping when a
  # worker fails; surface the first real error rather than a cryptic rbind
  bad <- !vapply(res, is.numeric, logical(1))
  if (any(bad)) stop("run_methods failed in ", sum(bad), " replicate(s): ",
                     as.character(res[[which(bad)[1]]]))
  do.call(rbind, res)
}

pvals_path <- function(cfg, n1, id)
  file.path(cfg$results_dir, sprintf("n1_%d", n1), sprintf("pvals_%04d.rds", id))

# The p-values depend on the data and on the Dunn adjustment.
pvals_key <- function(cfg) paste(data_key(cfg), cfg$dunn_p_adjust)

# Run the methods for the given scenario ids of one n1, saving one p-value
# file per scenario. Existing current files are skipped (restart after a crash).
# Data are loaded one scenario at a time to keep memory low.
run_methods_chunk <- function(cfg, n1, ids) {
  key <- pvals_key(cfg)
  dir.create(dirname(pvals_path(cfg, n1, 1)), showWarnings = FALSE, recursive = TRUE)
  for (id in ids) {
    f <- pvals_path(cfg, n1, id)
    if (file.exists(f) && identical(readRDS(f)$key, key)) next
    x <- load_scenario_data(cfg, n1, id)
    d <- x$design
    p <- run_scenario(x$data, c(d$n1, d$n2, d$n3), cfg)
    saveRDS(list(key = key, design = d, pvals = p), f)
  }
  invisible(NULL)
}

# Read all saved p-value files back into list(design, pvals), one element per
# n1, in the form extract_rejections() expects. Errors if any are missing.
load_results <- function(cfg) {
  key <- pvals_key(cfg)
  lapply(cfg$n1_grid, function(n1) {
    ids <- design_for_n1(cfg, n1)$scenario_id
    files <- vapply(ids, function(id) pvals_path(cfg, n1, id), "")
    missing <- !file.exists(files)
    if (any(missing)) stop(sum(missing), " of ", length(ids),
                           " p-value files missing for n1 = ", n1,
                           " (first: ", files[which(missing)[1]], ")")
    x <- lapply(files, readRDS)
    if (!all(vapply(x, function(z) identical(z$key, key), logical(1))))
      stop("Stale p-value files for n1 = ", n1, " (settings changed)")
    list(design = do.call(rbind, lapply(x, `[[`, "design")),
         pvals  = lapply(x, `[[`, "pvals"))
  })
}
