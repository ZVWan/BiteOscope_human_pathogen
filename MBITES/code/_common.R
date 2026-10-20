# Shared setup for the runs.
#
# The model code (configuration, landscape and the infection-dependent
# manipulation) sits in code/manuscript/ as a copy, so this folder runs without the
# repository beside it; check_model_snapshot() stops a run if those copies change.
# On top of it:
#
#   EIP_DAYS <- 11L   the extrinsic incubation period, used for the simulated
#                     biology and for the vectorial capacity readout, so that the
#                     two describe the same process.
#
# pin_model_eip() writes the package value explicitly instead of relying on the
# default, and check_units() fails loudly if model and readout ever disagree.

EIP_DAYS <- 11L                                  # the reported incubation period

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
.here <- if (exists("SCRIPT_DIR")) SCRIPT_DIR else .fb_dir()
RUN_ROOT <- normalizePath(file.path(.here, ".."))
MS_CODE <- normalizePath(file.path(RUN_ROOT, "code", "manuscript"), mustWork = FALSE)
if (!file.exists(file.path(MS_CODE, "_senescence-common.R")))
  stop("model code snapshot not found at ", MS_CODE)
rm(.here)

# ---- model code snapshot -----------------------------------------------------
# A copy of the model code in the repository, with the comments shortened and the
# incubation period set to 11 d. The hashes are the revision that produced the
# stored results: a changed copy stops the run instead of quietly scoring a
# different model. To move to a newer revision, copy both files in again and paste
# the hashes the message prints.

MANUSCRIPT_MD5 <- c(
  "_senescence-common.R"            = "02337212415a135841a3efd9ab812653",
  "manipulation-enhanced-probing.R" = "9d3b6e1db87d4cb4415955b62a8744bf"
)

check_model_snapshot <- function(dir = MS_CODE) {
  got <- unname(tools::md5sum(file.path(dir, names(MANUSCRIPT_MD5))))
  if (anyNA(got) || !identical(got, unname(MANUSCRIPT_MD5)))
    stop("the model code in ", dir, " is not the revision the stored results came from.\n",
         "  found ", paste(sprintf('"%s" = "%s"', names(MANUSCRIPT_MD5), got),
                          collapse = ",\n        "),
         "\n  if the change is intentional, paste those values into MANUSCRIPT_MD5 ",
         "in code/_common.R.")
  invisible(TRUE)
}

check_model_snapshot()

Sys.setenv(MBITES_COMMON_DIR = MS_CODE)
suppressMessages(source(file.path(MS_CODE, "_senescence-common.R")))
EIP_DAYS <- 11L                                  # set after the copy is sourced

# the manuscript settings, not the small defaults in the common file
nMosquitos <- 1000L; lambda <- 110; nHumans <- 200L; init_prev <- 0.1

# ---- EIP ---------------------------------------------------------------------

pin_model_eip <- function(eip_days = EIP_DAYS) {
  priv <- MBITES:::PathogenParameters$.__enclos_env__$private
  priv$mosquito_incubation <- as.integer(eip_days)
  invisible(NULL)
}

check_units <- function(eip_days = EIP_DAYS) {
  got <- MBITES:::PathogenParameters$get_mosquito_incubation()
  if (!identical(as.integer(got), as.integer(eip_days)) || !identical(as.integer(EIP_DAYS), as.integer(eip_days)))
    stop("model EIP and readout EIP differ: model ", got, " d, readout ", EIP_DAYS, " d")
  invisible(NULL)
}

pin_model_eip()
check_units()

# ---- runs --------------------------------------------------------------------

# as run_gj() in the copied model code, with the EIP pinned before the run and
# the VC read at that same value. init_age_max keeps the manuscript behaviour
# (initial human infections aged 1-200 d) unless a short infection is wanted.
run_sim <- function(seed, manip = list(), init_age_max = 200L) {
  pin_model_eip()
  MBITES:::Globals$hardreset(); configure_model(); trackBloodHost()
  options(mbites.manip.probe_multiplier=1.0, mbites.manip.feed_multiplier=1.0,
          mbites.manip.surviveprobe_multiplier=1.0, mbites.manip.feed_survival=1.0,
          mbites.manip.bloodbout_rate_multiplier=1.0,
          mbites.manip.bloodattempt_rate_multiplier=1.0, mbites.manip.bloodsearch_rate_multiplier=1.0)
  if (length(manip)) options(manip)
  reset_manip_telemetry(); set.seed(seed)
  Tile_Initialize(make_landscape())
  hpr <- rbinom(nHumans, 1L, init_prev)
  Human_NULL_Initialize(data.frame(tileID=rep(1,nHumans), siteID=rep(1:2,length.out=nHumans),
                                    feedingID=rep(1,nHumans), w=0.9, pr=hpr))
  MBITES_Initialize(data.frame(tileID=rep(1,nMosquitos), siteID=rep(1:2,length.out=nMosquitos), female=TRUE))
  hh <- MBITES:::Globals$get_tile(1)$get_humans(); anc <- SEI_Pathogen$new(); aid <- anc$get_id()
  for (kk in hh$ls()) if (isTRUE(hpr[as.integer(kk)] == 1L)) {
    p <- SEI_Pathogen$new(aid); p$initial_age(sample(seq_len(init_age_max), 1)); hh$get(kk)$add_pathogen(p) }
  d <- paste0(file.path(tempdir(), sprintf("mbites_%d_%s", seed, paste(unlist(manip), collapse="_"))), "/")
  unlink(d, recursive = TRUE); set_output(directory = d, runID = 1)
  simulation(tMax = duration_days, pretty = TRUE, cleanup = TRUE)
  vc <- vc_readout(d); unlink(d, recursive = TRUE)
  data.frame(eip = EIP_DAYS, gj_vc = vc,
             inf_probes = getOption("mbites.manip.telemetry")$probes_infectious)
}

# Garrett-Jones VC at the reported EIP, from the bite history written by the run
vc_readout <- function(dir) {
  f <- file.path(dir, "run1", "mosquito_F_1.json")
  m <- jsonlite::fromJSON(f, flatten = TRUE)
  nn <- which(sapply(m$id, is.null)); if (length(nn)) m <- m[-nn, ]
  h <- data.frame(id = 1:nHumans)
  invisible(capture.output(v <- Bionomics_vectorialCapacity(m, h, EIP = EIP_DAYS)))
  sum(sapply(v, function(x) x$VC))
}

# paired percentile bootstrap on the pooled fold, as in the manuscript scripts.
# Returns NA when no seed pair has a usable denominator (an all-zero baseline).
pooled_fold <- function(y, x, B = 2000, seed = 20260908) {
  ok <- is.finite(y) & is.finite(x) & x > 0
  y <- y[ok]; x <- x[ok]
  if (!length(x) || sum(x) <= 0) return(c(est = NA_real_, lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  n <- length(x)
  b <- replicate(B, { k <- sample.int(n, replace = TRUE); sum(y[k]) / sum(x[k]) })
  q <- quantile(b, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  c(est = sum(y) / sum(x), lo = q[1], hi = q[2])
}

# ---- storage and driver ------------------------------------------------------
#
# A results file grows as the runs finish: every finished run is appended, and the
# file gets its DONE line only once all the tasks in it are done. An interrupted
# job therefore keeps everything it has already produced, and running the same
# script again continues where it stopped instead of starting over.

DONE_MARK <- "DONE"

# the rows already in a results file, or NULL if there is nothing stored yet
read_stored <- function(csv) {
  if (!file.exists(csv)) return(NULL)
  d <- utils::read.csv(csv, stringsAsFactors = FALSE, fill = TRUE)
  eip <- suppressWarnings(as.integer(d$eip))
  d <- d[!is.na(eip), , drop = FALSE]        # drops the trailing DONE line
  if (!nrow(d)) return(NULL)
  d$eip <- as.integer(d$eip)
  d$gj_vc <- as.numeric(d$gj_vc)
  if (!is.null(d$inf_probes)) d$inf_probes <- as.integer(d$inf_probes)
  if (!is.null(d$seed)) d$seed <- as.integer(d$seed)
  if (!is.null(d$mult)) d$mult <- as.numeric(d$mult)
  if (!is.null(d$feed_survival)) d$feed_survival <- as.numeric(d$feed_survival)
  d
}

is_complete <- function(csv) {
  if (!file.exists(csv)) return(FALSE)
  x <- readLines(csv, warn = FALSE)
  length(x) > 0L && identical(trimws(x[length(x)]), DONE_MARK)
}

# what the analysis scripts use: read a finished file, or stop on a partial one
read_complete <- function(csv) {
  if (!file.exists(csv)) stop("missing ", csv, " - run the matching runner first.")
  if (!is_complete(csv))
    stop("incomplete ", csv, ": no DONE line, so the runner was interrupted. ",
         "Run it again to finish the missing runs.")
  d <- read_stored(csv)
  if (is.null(d)) stop("no rows in ", csv)
  ids <- do.call(paste, c(d[setdiff(names(d), c("gj_vc", "inf_probes"))], sep = "\r"))
  if (any(duplicated(ids)))
    stop(sprintf("%s holds %d duplicate row(s), so some runs were stored twice. ",
                 csv, sum(duplicated(ids))),
         "Run the matching runner again: it keeps one copy of each and rewrites the file.")
  d
}

append_rows <- function(csv, rows) {
  first <- !file.exists(csv)
  utils::write.table(rows, csv, sep = ",", row.names = FALSE, col.names = first,
                     quote = FALSE, append = !first)
}

# overwrite a store with the given rows (used to drop half-written runs)
write_rows <- function(csv, rows) {
  if (!nrow(rows)) { unlink(csv); return(invisible(NULL)) }
  utils::write.table(rows, csv, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

mark_done <- function(csv) cat(DONE_MARK, "\n", file = csv, append = TRUE)

# PSOCK cluster with the model code snapshot and the baseline probing rate loaded
# on every worker; file is sourced relative to code/ with SCRIPT_DIR set, so the
# DENV runs get code/denv-pathogen.R instead
worker_cluster <- function(nw, root, script_subdir, file = "_common.R") {
  cl <- parallel::makeCluster(nw, type = "PSOCK")
  w_root <- root; w_subdir <- script_subdir; w_file <- file
  parallel::clusterExport(cl, c("w_root", "w_subdir", "w_file"), envir = environment())
  parallel::clusterEvalQ(cl, {
    SCRIPT_DIR <- file.path(w_root, w_subdir)
    suppressMessages(source(file.path(w_root, "code", w_file)))
    BASE_PROBEH <- 0.4
    rm(w_root, w_subdir, w_file)
  })
  cl
}

# Run a task table on the cluster, appending each finished chunk of runs to csv.
# fun(i) returns the row of task i; the task columns are the key, so tasks whose
# rows are already stored are skipped. Returns everything stored in csv.
run_tasks <- function(cl, tasks, csv, fun, label = basename(csv),
                          exports = character(), chunk_factor = 4L) {
  # numeric keys are formatted element by element: a value has to survive its round
  # trip through the CSV, and formatting a whole vector pads to a common width
  as_key <- function(v) {
    if (is.numeric(v)) vapply(v, function(z) format(z, digits = 12, trim = TRUE), character(1))
    else as.character(v)
  }
  key_of <- function(d) do.call(paste, c(lapply(d[, names(tasks), drop = FALSE], as_key),
                                         sep = "\r"))
  total <- nrow(tasks)
  done <- read_stored(csv)
  if (!is.null(done)) {
    # a run is done when its row is in the file; one without a row (killed
    # mid-write) is dropped and done again, so that rows are never stored twice
    k <- key_of(done)
    # two processes appending to the same file leave every run twice; keep one copy
    dup <- duplicated(paste(k, done$eip))
    if (any(dup)) {
      cat(sprintf("[%s] %d duplicate row(s) in %s, keeping one copy of each\n",
                  label, sum(dup), basename(csv)))
      done <- done[!dup, , drop = FALSE]
      k <- k[!dup]
      write_rows(csv, done)
    }
    has_row <- tapply(done$eip, k, function(e) EIP_DAYS %in% e)
    incomplete <- names(has_row)[!has_row]
    if (length(incomplete)) {
      cat(sprintf("[%s] %d run(s) in %s have no row, they will be done again\n",
                  label, length(incomplete), basename(csv)))
      keep <- !(k %in% incomplete)
      write_rows(csv, done[keep, , drop = FALSE])
      done <- done[keep, , drop = FALSE]
    }
    tasks <- tasks[!(key_of(tasks) %in% key_of(done)), , drop = FALSE]
    if (total > nrow(tasks))
      cat(sprintf("[%s] %d of %d runs already in %s -> skipped\n",
                  label, total - nrow(tasks), total, basename(csv)))
  }
  if (!nrow(tasks)) {
    if (!is_complete(csv)) mark_done(csv)
    cat(sprintf("[%s] complete: %d runs\n", label, total))
    return(invisible(read_stored(csv)))
  }

  cat(sprintf("[%s] START | workers=%d | runs=%d | EIP %d d\n",
              label, length(cl), nrow(tasks), EIP_DAYS))
  if (length(exports)) parallel::clusterExport(cl, exports, envir = parent.frame())
  # the workers must see the tasks that are actually left to do, in the same order
  # as the indices handed out below
  parallel::clusterExport(cl, "tasks", envir = environment())
  chunk <- max(1L, chunk_factor * length(cl))
  idx <- seq_len(nrow(tasks))
  t0 <- Sys.time(); n_done <- 0L; n_bad <- 0L
  for (g in split(idx, ceiling(idx / chunk))) {
    rows <- do.call(rbind, parallel::parLapply(cl, g, fun))
    append_rows(csv, rows)
    n_done <- n_done + length(g)
    n_bad <- n_bad + sum(is.na(rows$gj_vc))
    cat(sprintf("  [%s] %d/%d runs stored | %.1f min\n", label, n_done, nrow(tasks),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  mark_done(csv)
  if (n_bad > 0) cat(sprintf("[%s] WARNING: %d run(s) errored (gj_vc = NA)\n", label, n_bad))
  cat(sprintf("[%s] DONE in %.1f min -> %s\n", label,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), csv))
  invisible(read_stored(csv))
}
