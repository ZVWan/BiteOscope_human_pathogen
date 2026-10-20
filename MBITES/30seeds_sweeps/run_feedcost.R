# Fig. 4d trade-off (S6): the engorgement boost with and without a
# per-meal survival cost, in one run.
#
#   boost arm     scenario "feedcost"          feedH x1.92 AND feed_survival = fs
#                 scenario "baseline"          no manipulation (fold denominator)
#   no-boost arm  scenario "feedcost_noboost"  feed_survival = fs only
#
# Usage: Rscript run_feedcost.R [n_workers]
#
# FS_GRID runs from 1.0 (no cost) down to 0.70, i.e. 0-30% mortality per meal;
# fs = 1.0 in the boost arm is the plain S5 engorgement boost. Each row is one run
# at the reported EIP. Runs are appended to the CSV as they finish, so an
# interrupted job keeps what it has produced and running the script again
# continues. Output: results/gj-feedcost-30seeds.csv.

N_WORKERS <- 8
FEED_OR   <- 1.92
FS_GRID   <- seq(1.0, 0.70, by = -0.05)

.fb_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (length(p) && nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
SCRIPT_DIR <- .fb_dir()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
OUT <- file.path(SCRIPT_DIR, "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CODE_DIR <- file.path(ROOT, "code")

args <- commandArgs(trailingOnly = TRUE)
nw_arg <- if (length(args) >= 1) args[1] else NA_character_
nw <- if (!is.na(nw_arg)) suppressWarnings(as.integer(nw_arg)) else N_WORKERS
if (is.na(nw) || nw < 1) nw <- min(parallel::detectCores() - 1, 8)

suppressMessages(source(file.path(CODE_DIR, "_common.R")))
BASE_PROBEH <- 0.4

SEEDS30 <- seq(101L, by = 101L, length.out = 30L)
tasks <- rbind(
  data.frame(arm = "boost",    scenario = "baseline",          seed = SEEDS30, feed_survival = 1.0),
  expand.grid(arm = "boost",    scenario = "feedcost",          seed = SEEDS30,
              feed_survival = FS_GRID, stringsAsFactors = FALSE),
  expand.grid(arm = "noboost",  scenario = "feedcost_noboost",  seed = SEEDS30,
              feed_survival = FS_GRID, stringsAsFactors = FALSE))
tasks <- tasks[order(tasks$arm, tasks$seed, tasks$feed_survival), ]
csv <- file.path(OUT, "gj-feedcost-30seeds.csv")

run_one <- function(i) {
  nm <- tasks$scenario[i]; s <- tasks$seed[i]; fs <- tasks$feed_survival[i]
  arm <- tasks$arm[i]
  manip <- if (nm == "baseline") list() else
    if (nm == "feedcost") list(mbites.manip.feed_multiplier = FEED_OR,
                               mbites.manip.feed_survival = fs) else
      list(mbites.manip.feed_survival = fs)
  tryCatch({
    cbind(arm = arm, scenario = nm, feed_survival = fs, seed = s, run_sim(s, manip))
  }, error = function(e) {
    data.frame(arm = arm, scenario = nm, feed_survival = fs, seed = s,
               eip = EIP_DAYS, gj_vc = NA_real_, inf_probes = NA_integer_)
  })
}

cl <- worker_cluster(nw, ROOT, "30seeds_sweeps")
out <- run_tasks(cl, tasks, csv, run_one, label = "S6 trade-off", exports = "FEED_OR")
parallel::stopCluster(cl)

# pooled fold per fs level: boost against "baseline",
# no-boost against its own fs = 1.0 rows
z <- out[out$eip == EIP_DAYS & !is.na(out$gj_vc), ]
base_boost <- z[z$scenario == "baseline", c("seed", "gj_vc")]; names(base_boost)[2] <- "b"
base_noboost <- z[z$scenario == "feedcost_noboost" & abs(z$feed_survival - 1) < 1e-9,
                  c("seed", "gj_vc")]; names(base_noboost)[2] <- "b"
cat(sprintf("\npooled fold vs no-manipulation baseline at EIP %d d (per-meal death = 1 - fs):\n", EIP_DAYS))
for (arm in c("boost", "noboost")) {
  nm <- if (arm == "boost") "feedcost" else "feedcost_noboost"
  base <- if (arm == "boost") base_boost else base_noboost
  for (fs in sort(unique(z$feed_survival[z$scenario == nm]), decreasing = TRUE)) {
    a <- z[z$scenario == nm & abs(z$feed_survival - fs) < 1e-9, c("seed", "gj_vc")]
    m <- merge(a, base, by = "seed")
    ci <- pooled_fold(m$gj_vc, m$b)
    cat(sprintf("  %-8s fs %.2f (%.0f%% death)  %.3f [%.3f, %.3f]\n",
                arm, fs, 100 * (1 - fs), ci[["est"]], ci[["lo"]], ci[["hi"]]))
  }
}
