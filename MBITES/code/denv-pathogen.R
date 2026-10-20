# Dengue pathogen settings for the denv_sweeps scripts, on top of the shared setup.
#
# MBITES fixes five pathogen values in its own source. Four of them are reachable
# at run time and are set here; the incubation period stays at 11 d (see _common.R), so only
# the transmission and human-infection side changes.
#
#   setting                value  package  source
#   b (mosquito -> human)   0.55     0.55  Andraud et al. 2012 PLoS ONE review:
#                                          dengue models use 0.1-0.75.
#   c (human -> mosquito)   0.50     0.15  Andraud et al. 2012: 0.5-1.
#   human incubation         6 d      6 d  Chan & Johansson 2012 PLoS ONE: mean
#                                          5.9 d (95% 3-10 d).
#   human infectious        5 d    194 d  WHO: viraemic about 4-5 d. Nguyen et al.
#                                          2013 PNAS: infectious 4-5 d after onset.
#
# The package writes the human clearance (200 d) directly into the body of
# SEI_Pathogen$oneDay_human(), where no getter or setter reaches it, so that
# method is rebound from the copy below. Sourcing this file changes the pathogen
# settings for the rest of the session, so keep it out of the manuscript runs.

DENV_B <- 0.55
DENV_C <- 0.50
DENV_HUMAN_INCUBATION_DAYS <- 6L
DENV_HUMAN_INFECTIOUS_DAYS <- 5L
DENV_HUMAN_CLEARANCE_DAYS <- DENV_HUMAN_INCUBATION_DAYS + DENV_HUMAN_INFECTIOUS_DAYS

# requires the caller to have set SCRIPT_DIR (every runner and plot script in this
# folder does), which is how the shared code finds the manuscript folder
if (!exists("SCRIPT_DIR")) stop("SCRIPT_DIR is not set: source this from a runner or plot script in this folder")
source(file.path(SCRIPT_DIR, "..", "code", "_common.R"))

apply_denv_settings <- function() {
  priv <- MBITES:::PathogenParameters$.__enclos_env__$private
  priv$b <- DENV_B
  priv$c <- DENV_C
  priv$human_incubation <- DENV_HUMAN_INCUBATION_DAYS
  pin_model_eip(EIP_DAYS)
  invisible(NULL)
}

check_denv_settings <- function() {
  pp <- MBITES:::PathogenParameters
  ok <- pp$get_b() == DENV_B && pp$get_c() == DENV_C &&
    pp$get_human_incubation() == DENV_HUMAN_INCUBATION_DAYS &&
    pp$get_mosquito_incubation() == EIP_DAYS
  if (!ok)
    stop("DENV settings did not stick (b=", pp$get_b(), " c=", pp$get_c(),
         " IIP=", pp$get_human_incubation(), " EIP=", pp$get_mosquito_incubation(), ")")
  invisible(NULL)
}

# the package method with the 200 d clearance replaced by DENV_HUMAN_CLEARANCE_DAYS
denv_oneDay_human <- function(human) {
  if (is.null(private$time_created)) private$time_created <- MBITES:::Globals$get_tNow()
  private$incubating <- private$incubating + 1L
  if (private$incubating > DENV_HUMAN_CLEARANCE_DAYS) {
    human$remove_pathogen(self)
    self$exit()
  } else if (!private$infectious && private$incubating >= private$incubation_h) {
    private$infectious <- TRUE
    private$time_became_infectious <- MBITES:::Globals$get_tNow()
  }
}
SEI_Pathogen$set(which = "public", name = "oneDay_human",
                 value = denv_oneDay_human, overwrite = TRUE)

apply_denv_settings()
check_denv_settings()
cat(sprintf("[denv] b %.2f, c %.2f, human incubation %d d, infectious %d d (clearance day %d), EIP %d d\n",
            DENV_B, DENV_C, DENV_HUMAN_INCUBATION_DAYS, DENV_HUMAN_INFECTIOUS_DAYS,
            DENV_HUMAN_CLEARANCE_DAYS, EIP_DAYS))

# run_sim() with the DENV settings applied and the initial human infections aged
# inside the 11-day infection window instead of the 1-200 d default
run_denv <- function(seed, manip = list()) {
  apply_denv_settings()
  run_sim(seed, manip, init_age_max = DENV_HUMAN_CLEARANCE_DAYS)
}
