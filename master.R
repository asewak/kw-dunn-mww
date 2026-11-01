# master.R -- global parameters; reads/generates datasets, runs methods, plots.
#
# Usage (from the repo root):
#   Rscript master.R                       # everything, locally
#   Rscript master.R mode=run n1=100 chunk=3 nchunks=10
#                                          # cluster job: scenarios in chunk 3 of 10, n1 = 100
#   Rscript master.R mode=combine          # after all jobs: rejection rates, CSV, plots
# Arguments (all optional): mode = all | run | combine; n1 = one n1 value
# (default all in n1_grid); chunk / nchunks = split an n1's scenarios across jobs.
# Existing current data / p-value files are skipped, so a killed job can be
# simply resubmitted.

cfg <- list(
  reps        = 5,                         # replications per scenario (small for testing)
  n1_grid     = c(5, 15, 30, 50, 100),     # n1 = n2
  ratio_grid  = c(0.5, 1, 2, 4),           # n3 / n1
  mu3_grid    = c(0, 0.5, 1, 2, 3),        # location of group 3 (location, skewed)
  sigma3_grid = c(0.1, 0.25, 0.5, 1, 2, 4),# scale of group 3 (scale family)
  delta_grid  = c(0, 0.3, 0.5, 0.8),       # shift of group 2 vs 1 (0 = null)
  alpha       = 0.05,
  dunn_p_adjust = "none",                  # adjustment across the 3 Dunn pairs
  # cores for mclapply: SLURM allocation if present, else 4
  mc_cores    = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4")),
  parallel_min_reps = 200,                 # below this many reps, run serially (fork overhead)
  seed        = 2024,
  data_dir    = "data",
  results_dir = "results",
  output_dir  = "output"
)

source("R/dgps.R")
source("R/methods.R")
source("R/outputs.R")

# ---- command-line arguments ------------------------------------------------
args <- strsplit(commandArgs(trailingOnly = TRUE), "=", fixed = TRUE)
arg <- function(name, default) {
  hit <- Filter(function(a) a[1] == name, args)
  if (length(hit)) hit[[1]][2] else default
}
mode    <- arg("mode", "all")
n1_run  <- if (is.na(n1 <- suppressWarnings(as.numeric(arg("n1", NA))))) cfg$n1_grid else n1
chunk   <- as.integer(arg("chunk", 1))
nchunks <- as.integer(arg("nchunks", 1))
stopifnot(mode %in% c("all", "run", "combine"), chunk >= 1, chunk <= nchunks)

# ---- 1 & 2. data + methods (per scenario files; restartable) ----------------
if (mode %in% c("all", "run")) {
  for (n1 in n1_run) {
    ids <- design_for_n1(cfg, n1)$scenario_id
    ids <- ids[(seq_along(ids) - 1) %% nchunks == chunk - 1]
    message("n1 = ", n1, ": ", length(ids), " scenarios (chunk ", chunk, "/", nchunks, ")")
    save_scenario_data(cfg, n1, ids)     # DGPs: generated once, then reused
    run_methods_chunk(cfg, n1, ids)      # methods: p-values saved per scenario
  }
}

# ---- 3. outputs ------------------------------------------------------------
if (mode %in% c("all", "combine")) {
  results <- load_results(cfg)
  saveRDS(list(cfg = cfg, results = results),
          file.path(cfg$results_dir, "all_results.rds"))
  rej <- extract_rejections(results, cfg)
  saveRDS(rej, file.path(cfg$results_dir, "rejection_rates.rds"))
  save_rejection_csv(rej, cfg)
  save_provenance(cfg)
  save_all_plots(rej, cfg)
  message("Done. Plots and CSV in ", cfg$output_dir)
}
