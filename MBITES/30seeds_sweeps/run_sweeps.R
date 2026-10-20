# The three Garrett-Jones sweeps at 30 seeds:
#   probing      probeH             -> results/gj-probing-30seeds.csv         (Fig. 4b)
#   landing rate bloodattempt rate  -> results/gj-rateb-probeH04-30seeds.csv  (Fig. 4a)
#   engorgement  feedH              -> results/gj-feed-30seeds.csv            (Fig. 4c)
#
# Usage: Rscript run_sweeps.R [n_workers]
#
# All three run at the common probing baseline probeH 0.4. Each row is one run at
# the reported EIP. Runs are appended to the CSV as they finish, so an interrupted
# job keeps what it has produced and running the script again continues.

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

SEEDS30 <- seq(101L, by = 101L, length.out = 30L)
grid <- seq(1, 2.5, by = 0.25)

sweeps <- list(
  probing = list(label = "probing",      knob = "probe_multiplier",
                 csv = "gj-probing-30seeds.csv"),
  rateb04 = list(label = "landing rate", knob = "bloodattempt_rate_multiplier",
                 csv = "gj-rateb-probeH04-30seeds.csv"),
  feed    = list(label = "engorgement",  knob = "feed_multiplier",
                 csv = "gj-feed-30seeds.csv")
)

run_one <- function(i) {
  s <- tasks$seed[i]; m <- tasks$mult[i]
  tryCatch({
    manip <- list(); manip[[paste0("mbites.manip.", knob)]] <- m
    cbind(seed = s, mult = m, run_sim(s, manip))
  }, error = function(e) {
    data.frame(seed = s, mult = m, eip = EIP_DAYS, gj_vc = NA_real_, inf_probes = NA_integer_)
  })
}

cl <- worker_cluster(nw, ROOT, "30seeds_sweeps")
for (id in names(sweeps)) {
  def <- sweeps[[id]]
  knob <- def$knob
  tasks <- expand.grid(seed = SEEDS30, mult = grid, stringsAsFactors = FALSE)
  tasks <- tasks[order(tasks$seed, tasks$mult), ]
  run_tasks(cl, tasks, file.path(OUT, def$csv), run_one, label = def$label,
                exports = "knob")
}
parallel::stopCluster(cl)
