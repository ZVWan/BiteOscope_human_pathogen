# Full Fig. 4 rerun with the per-infected-mosquito readout (metric-counter.R), same
# seeds and grids as the stored runs:
#   malaria settings  probing / rate_b / engorgement sweeps (1-2.5 x), S6 feed cost
#   dengue settings   probing / rate_b sweeps
# The exact points are already in results/metric-points-30seeds.csv.
#
# Usage: Rscript run_metric_all.R [n_workers]
# Stopping: create a file named STOP in this folder (touch STOP); the run halts after
# the current batch (at most ~2 min). Starting the script again skips every run
# already stored in results/metric-sweeps-30seeds.csv and removes STOP.
args <- commandArgs(trailingOnly = TRUE)
nw <- if (length(args)) as.integer(args[1]) else 12L

# Paths are resolved from this script, so the package runs from any checkout.
# MBITES_MODEL_DIR can point at another copy of the model tree.
SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (length(p) && nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
})
HERE <- normalizePath(SCRIPT_DIR)
EN   <- normalizePath(Sys.getenv("MBITES_MODEL_DIR", file.path(HERE, "..", "MBITES")))
STOP <- file.path(HERE, "STOP"); unlink(STOP)
csv  <- file.path(HERE, "results", "metric-sweeps-30seeds.csv")
SEEDS30 <- seq(101L, by = 101L, length.out = 30L)
GRID    <- seq(1, 2.5, by = 0.25)
FS_GRID <- seq(1.0, 0.70, by = -0.05)
KNOB <- c(probing = "mbites.manip.probe_multiplier",
          rateb   = "mbites.manip.bloodattempt_rate_multiplier",
          feed    = "mbites.manip.feed_multiplier")

sweep_tasks <- function(set, studies)
  do.call(rbind, lapply(studies, function(st)
    data.frame(set = set, study = st, scenario = st, x = rep(GRID, each = 30),
               seed = SEEDS30, stringsAsFactors = FALSE)))
fc <- rbind(
  data.frame(scenario = "baseline", x = 1.0, seed = SEEDS30),
  expand.grid(scenario = "feedcost", x = FS_GRID, seed = SEEDS30, stringsAsFactors = FALSE),
  expand.grid(scenario = "feedcost_noboost", x = FS_GRID, seed = SEEDS30, stringsAsFactors = FALSE))
fc <- cbind(set = "malaria", study = "feedcost", fc, stringsAsFactors = FALSE)

sets <- list(
  malaria = list(file = "_common.R", runner = "run_sim",
                 tasks = rbind(sweep_tasks("malaria", c("rateb", "probing", "feed")), fc)),
  dengue  = list(file = "denv-pathogen.R", runner = "run_denv",
                 tasks = sweep_tasks("dengue", c("rateb", "probing"))))

manip_of <- function(study, scenario, x) {
  if (study != "feedcost") return(setNames(list(x), KNOB[[study]]))
  switch(scenario, baseline = list(),
         feedcost = list(mbites.manip.feed_multiplier = 1.92, mbites.manip.feed_survival = x),
         feedcost_noboost = list(mbites.manip.feed_survival = x))
}
key <- function(d) paste(d$set, d$study, d$scenario, sprintf("%.2f", d$x), d$seed)

for (set in names(sets)) {
  def <- sets[[set]]; tasks <- def$tasks
  done <- if (file.exists(csv)) read.csv(csv, stringsAsFactors = FALSE) else NULL
  if (!is.null(done)) tasks <- tasks[!key(tasks) %in% key(done), ]
  # interleave seeds so every sweep fills in evenly if the run is stopped early
  tasks <- tasks[order(match(tasks$seed, SEEDS30), tasks$study, tasks$scenario, tasks$x), ]
  cat(sprintf("[%s] %d runs to do\n", set, nrow(tasks)))
  if (!nrow(tasks)) next
  if (file.exists(STOP)) break
  cl <- parallel::makeCluster(nw, type = "PSOCK")
  parallel::clusterExport(cl, c("EN", "HERE", "def", "set", "tasks", "manip_of", "KNOB"))
  parallel::clusterEvalQ(cl, {
    SCRIPT_DIR <- EN
    suppressMessages(source(file.path(EN, "code", def$file)))
    BASE_PROBEH <- 0.4
    source(file.path(HERE, "metric-counter.R"))
  })
  run_one <- function(i) {
    t <- tasks[i, ]
    tryCatch({
      reset_metric()
      r <- get(def$runner)(t$seed, manip_of(t$study, t$scenario, t$x))
      cbind(t, r, metric_summary())
    }, error = function(e) {
      message("run failed: ", paste(t, collapse = " "), ": ", conditionMessage(e)); NULL })
  }
  t0 <- Sys.time(); idx <- seq_len(nrow(tasks))
  for (g in split(idx, ceiling(idx / nw))) {
    if (file.exists(STOP)) { cat("STOP file found, halting; run the script again to resume\n"); break }
    rows <- do.call(rbind, parallel::parLapply(cl, g, run_one))
    if (!is.null(rows)) write.table(rows, csv, sep = ",", row.names = FALSE,
                                    col.names = !file.exists(csv), append = file.exists(csv))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] %d/%d | %.0f min | ~%.0f min left in this set\n", set, max(g), nrow(tasks),
                el, el / max(g) * (nrow(tasks) - max(g))))
  }
  parallel::stopCluster(cl)
}
cat(if (file.exists(STOP)) "STOPPED\n" else "DONE\n")
