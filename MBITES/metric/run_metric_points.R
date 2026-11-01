# Test of the per-infected-mosquito readout on the Fig. 4 exact points:
# malaria settings (baseline04, S2_rate04, S3, S4, S5) and dengue settings
# (baseline, S2_rate, S4), 30 seeds each, same seeds and settings as the stored runs.
# Usage: Rscript run_metric_points.R [n_workers]
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
SEEDS30 <- seq(101L, by = 101L, length.out = 30L)

sets <- list(
  malaria = list(file = "_common.R", runner = "run_sim",
    scen = list(baseline04 = list(),
                S2_rate04  = list(mbites.manip.bloodattempt_rate_multiplier = 1.6),
                S3_probe   = list(mbites.manip.probe_multiplier = 1.29),
                S4_probe   = list(mbites.manip.probe_multiplier = 1.52),
                S5_feed    = list(mbites.manip.feed_multiplier  = 1.92))),
  dengue = list(file = "denv-pathogen.R", runner = "run_denv",
    scen = list(baseline = list(),
                S2_rate  = list(mbites.manip.bloodattempt_rate_multiplier = 1.6),
                S4_probe = list(mbites.manip.probe_multiplier = 1.52)))
)

csv <- file.path(HERE, "results", "metric-points-30seeds.csv")
done <- if (file.exists(csv)) read.csv(csv) else NULL

for (set in names(sets)) {
  def <- sets[[set]]
  tasks <- expand.grid(seed = SEEDS30, scenario = names(def$scen), stringsAsFactors = FALSE)
  if (!is.null(done)) tasks <- tasks[!paste(set, tasks$scenario, tasks$seed) %in%
                                       paste(done$set, done$scenario, done$seed), ]
  if (!nrow(tasks)) next
  cl <- parallel::makeCluster(nw, type = "PSOCK")
  parallel::clusterExport(cl, c("EN", "HERE", "def", "set", "tasks"))
  parallel::clusterEvalQ(cl, {
    SCRIPT_DIR <- EN
    suppressMessages(source(file.path(EN, "code", def$file)))
    BASE_PROBEH <- 0.4
    source(file.path(HERE, "metric-counter.R"))
  })
  run_one <- function(i) {
    s <- tasks$seed[i]; nm <- tasks$scenario[i]
    tryCatch({
      reset_metric()
      r <- get(def$runner)(s, def$scen[[nm]])
      cbind(set = set, scenario = nm, seed = s, r, metric_summary())
    }, error = function(e) {
      message("run failed: ", set, " ", nm, " ", s, ": ", conditionMessage(e)); NULL })
  }
  t0 <- Sys.time()
  idx <- seq_len(nrow(tasks))
  for (g in split(idx, ceiling(idx / (2 * nw)))) {
    rows <- do.call(rbind, parallel::parLapply(cl, g, run_one))
    if (!is.null(rows)) write.table(rows, csv, sep = ",", row.names = FALSE,
                                    col.names = !file.exists(csv), append = file.exists(csv))
    cat(sprintf("[%s] %d/%d | %.1f min\n", set, max(g), nrow(tasks),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  parallel::stopCluster(cl)
}
cat("DONE\n")
