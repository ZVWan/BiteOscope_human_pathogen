# Shared setup for the model runs.
#
# The model code (configuration, landscape and the infection-dependent
# manipulation) sits in code/manuscript/ as a copy, so this folder runs without the
# repository beside it; check_model_snapshot() stops a run if those copies change.
#
# On top of it:
#
#   EIP_DAYS <- 11L        the extrinsic incubation period used by the simulation
#   run_sim(seed, manip)   one run, returns the telemetry counts
#
# pin_model_eip() writes the package value explicitly instead of relying on the
# default, and check_units() fails loudly if the model value drifts from EIP_DAYS.

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
# SCRIPT_DIR may be the model root (the folder that holds code/) or a folder inside it
RUN_ROOT <- if (file.exists(file.path(.here, "code", "manuscript"))) .here else
  normalizePath(file.path(.here, ".."))
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
  "_senescence-common.R"            = "a64c83e6bfa4de404625b41294876d06",
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

# The EIP is pinned before the run. init_age_max keeps the manuscript behaviour
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
  unlink(d, recursive = TRUE)
  data.frame(eip = EIP_DAYS,
             inf_probes = getOption("mbites.manip.telemetry")$probes_infectious)
}
