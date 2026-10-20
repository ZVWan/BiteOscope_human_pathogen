# The exact scenario points behind the Fig. 4 red dots, 30 seeds each:
#   baseline04   no manipulation, the shared denominator
#   S2_rate04    rate_b x1.6     DENV landing rate
#   S3_probe     probeH x1.29    NF135 probing
#   S4_probe     probeH x1.52    DENV probing
#   S5_feed      feedH x1.92     NF54 engorgement
#
# Usage: Rscript run_exact_points.R [n_workers]
#
# All five groups run at the probing baseline probeH 0.4. Each row is one run at
# the reported EIP. Runs are appended to the CSV as they finish, so an interrupted
# job keeps what it has produced and running the script again continues.
# Output: results/gj-exact-points-30seeds.csv.

N_WORKERS <- 8

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

scen <- list(
  baseline04 = list(),
  S2_rate04  = list(mbites.manip.bloodattempt_rate_multiplier = 1.6),
  S3_probe   = list(mbites.manip.probe_multiplier = 1.29),
  S4_probe   = list(mbites.manip.probe_multiplier = 1.52),
  S5_feed    = list(mbites.manip.feed_multiplier  = 1.92)
)
SEEDS30 <- seq(101L, by = 101L, length.out = 30L)

tasks <- do.call(rbind, lapply(names(scen), function(nm)
  data.frame(scenario = nm, seed = SEEDS30, stringsAsFactors = FALSE)))
tasks <- tasks[order(tasks$scenario, tasks$seed), ]
csv <- file.path(OUT, "gj-exact-points-30seeds.csv")

run_one <- function(i) {
  nm <- tasks$scenario[i]; s <- tasks$seed[i]
  tryCatch({
    cbind(scenario = nm, seed = s, run_sim(s, scen[[nm]]))
  }, error = function(e) {
    data.frame(scenario = nm, seed = s, eip = EIP_DAYS, gj_vc = NA_real_, inf_probes = NA_integer_)
  })
}

cl <- worker_cluster(nw, ROOT, "30seeds_sweeps")
out <- run_tasks(cl, tasks, csv, run_one, label = "exact points", exports = "scen")
parallel::stopCluster(cl)

# pooled fold per scenario, against baseline04 on the same seeds
z <- out[out$eip == EIP_DAYS & !is.na(out$gj_vc), ]
base <- z[z$scenario == "baseline04", c("seed", "gj_vc")]
names(base)[2] <- "b"
cat(sprintf("\npooled fold vs baseline04 at EIP %d d:\n", EIP_DAYS))
for (nm in names(scen)[-1]) {
  a <- z[z$scenario == nm, c("seed", "gj_vc")]
  m <- merge(a, base, by = "seed")
  ci <- pooled_fold(m$gj_vc, m$b)
  cat(sprintf("  %-11s %.3f [%.3f, %.3f] (n = %d)\n", nm, ci[["est"]], ci[["lo"]], ci[["hi"]], nrow(m)))
}
