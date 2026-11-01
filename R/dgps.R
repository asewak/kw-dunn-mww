# dgps.R -- data-generating processes. Generates and saves datasets (one file
# per scenario) so they are not regenerated on every run.
#
# Groups: 1 ~ N(0,1); 2 ~ N(delta,1); 3 depends on `family`:
#   location : N(mu3, 1)                      (sigma3 = 1)
#   scale    : N(0, sigma3^2)                 (mu3 = 0)
#   skewed   : mu3 + sigma3 * (Exp(1) - 1)    (standardised exponential, skew = 2)
# delta = 0 with mu3 = 0, sigma3 = 1 in the normal families is the global null;
# delta = 0 otherwise is a partial null (groups 1 & 2 equal, 3 differs);
# delta > 0 is the power setting.
# Sample sizes: n1 = n2, n3 = ceiling(ratio * n1).
#
# A scenario is identified by (n1, scenario_id). scenario_id is the row of the
# design table, which depends only on the grids in cfg. Each scenario has its
# own seed, so scenarios can be generated independently (e.g. in separate
# cluster jobs) and reproducibly.

make_design <- function(cfg) {
  with(cfg, {
    loc <- expand.grid(family = "location", mu3 = mu3_grid, sigma3 = 1,
                       stringsAsFactors = FALSE)
    scl <- expand.grid(family = "scale", mu3 = 0, sigma3 = sigma3_grid,
                       stringsAsFactors = FALSE)
    skw <- expand.grid(family = "skewed", mu3 = mu3_grid, sigma3 = 1,
                       stringsAsFactors = FALSE)
    g3 <- rbind(loc, scl, skw)
    d <- merge(g3, expand.grid(delta = delta_grid, ratio = ratio_grid))
    d <- d[order(d$family, d$delta, d$ratio, d$mu3, d$sigma3), ]
    d$scenario_id <- seq_len(nrow(d))
    rownames(d) <- NULL
    d
  })
}

# Design table for one n1 (adds the sample sizes).
design_for_n1 <- function(cfg, n1) {
  d <- make_design(cfg)
  d$n1 <- n1
  d$n2 <- n1
  d$n3 <- ceiling(d$ratio * n1)
  d
}

gen_group3 <- function(n, family, mu3, sigma3) {
  switch(family,
    location = rnorm(n, mu3, 1),
    scale    = rnorm(n, 0, sigma3),
    skewed   = mu3 + sigma3 * (rexp(n) - 1),
    stop("unknown family: ", family))
}

# One replicate: the stacked response vector (group labels are implied by
# n1, n2, n3 stored in the design table).
gen_one <- function(n1, n2, n3, delta, family, mu3, sigma3) {
  c(rnorm(n1), rnorm(n2, delta), gen_group3(n3, family, mu3, sigma3))
}

# Everything that determines the data. Saved with each file; a mismatch on
# load means the file is stale and it is regenerated.
data_key <- function(cfg) {
  paste(deparse(cfg[c("reps", "n1_grid", "ratio_grid", "mu3_grid",
                      "sigma3_grid", "delta_grid", "seed")]), collapse = "")
}

scenario_seed <- function(cfg, n1, id) cfg$seed + as.integer(n1) * 100000L + id

data_path <- function(cfg, n1, id)
  file.path(cfg$data_dir, sprintf("n1_%d", n1), sprintf("scen_%04d.rds", id))

# Generate the `reps` datasets of one design row.
generate_scenario <- function(cfg, r) {
  set.seed(scenario_seed(cfg, r$n1, r$scenario_id))
  replicate(cfg$reps,
            gen_one(r$n1, r$n2, r$n3, r$delta, r$family, r$mu3, r$sigma3),
            simplify = FALSE)
}

# Generate and save the data for the given scenario ids of one n1 (skips files
# that exist and are current). Returns invisibly.
save_scenario_data <- function(cfg, n1, ids, overwrite = FALSE) {
  d <- design_for_n1(cfg, n1)
  key <- data_key(cfg)
  dir.create(dirname(data_path(cfg, n1, 1)), showWarnings = FALSE, recursive = TRUE)
  for (id in ids) {
    f <- data_path(cfg, n1, id)
    if (!overwrite && file.exists(f) && identical(readRDS(f)$key, key)) next
    saveRDS(list(key = key, design = d[id, ],
                 data = generate_scenario(cfg, d[id, ])), f)
  }
  invisible(NULL)
}

load_scenario_data <- function(cfg, n1, id) {
  x <- readRDS(data_path(cfg, n1, id))
  if (!identical(x$key, data_key(cfg)))
    stop("Stale data file ", data_path(cfg, n1, id),
         " (settings changed); regenerate it")
  x
}
