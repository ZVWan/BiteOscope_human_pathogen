# Model settings and landscape shared by the runs.
#
# Sourced from code/_common.R, which sets MBITES_COMMON_DIR so that the
# infection-dependent manipulation layer is found from this copy.

suppressMessages({ library(MBITES); library(futile.logger) })
MBITES::local_logging("fatal", tempfile(fileext = ".log"))


# Load the manipulation framework; callers set MBITES_COMMON_DIR.
.mbites_common_dir <- Sys.getenv("MBITES_COMMON_DIR", unset = "")
if (nzchar(.mbites_common_dir)) {
  source(file.path(.mbites_common_dir, "manipulation-enhanced-probing.R"))
} else {
  source(file.path("scripts", "manipulation-enhanced-probing.R"))
}
rm(.mbites_common_dir)

# Shared settings; sweeps override these.
nMosquitos <- 150L; nHumans <- 200L; init_prev <- 0.5; lambda <- 8
duration_days <- 100
SNS_A <- 0.25; SNS_B <- 60          # age-based senescence; onset ~16 d
BASE_PROBEH <- 0.4                  # headroom for a 2.5x boost (0.4*2.5 = 1.0)
BASE_FEEDH  <- 0.5

make_landscape <- function() {
  L <- vector("list", 2)
  for (i in 1:2) { L[[i]]$id <- i; L[[i]]$xy <- c(1,1); L[[i]]$type <- 1L; L[[i]]$tileID <- 1L
    L[[i]]$move <- 0.5; L[[i]]$move_id <- (1:2)[-i]; L[[i]]$haz <- 0.005
    L[[i]]$feed <- list(list(w=1, enterP=1, zoo_id=-1, zoo_w=0.001))
    L[[i]]$aqua <- list(list(w=1, lambda=lambda)) }
  L
}

configure_model <- function() {
  MBITES::MBITES_Setup_Timing(timing_model=2, rate_b=1/(3/24),tmin_b=0, rate_bs=1/(6/24),tmin_bs=0,
    rate_o=1/(3/24),tmin_o=0, rate_os=1/(6/24),tmin_os=0, ppr_model=2, rate_ppr=1/(18/24), tmin_ppr=1e-2)
  MBITES:::Globals$set_SETUP("timing")
  Mosquito_Female$set("public","Overfeeding",MBITES:::mbites_OverfeedingNull,overwrite=TRUE); MBITES:::Globals$set_SETUP("bloodmeal")
  Mosquito_Female$set("public","rBatchSize",MBITES:::mbites_rBatchSizeBms,overwrite=TRUE)
  Mosquito_Female$set("public","Oogenesis",MBITES:::mbites_oogenesisMBDETES,overwrite=TRUE)
  Mosquito_Female$set("public","checkRefeed",MBITES:::mbites_checkRefeedMBDETES,overwrite=TRUE)
  Mosquito_Female$set("public","pReFeed",MBITES:::mbites_pReFeed_bm,overwrite=TRUE); MBITES:::Globals$set_SETUP("oogenesis")
  Mosquito_Female$set("public","queueSugarBout",MBITES:::mbites_queueSugarBout_null,overwrite=TRUE); MBITES:::Globals$set_SETUP("energetics")
  Mosquito_Female$set("public","layEggs",MBITES:::mbites_layEggs_Emerge,overwrite=TRUE); MBITES:::Globals$set_SETUP("oviposition")
  Mosquito$set("public","WingTattering",MBITES:::mbites_WingTattering_null,overwrite=TRUE)
  Mosquito$set("public","Senescence",MBITES:::mbites_Senescence,overwrite=TRUE)   # senescence on
  MBITES:::Globals$set_SETUP("survival")
  PATHOGEN_Setup("SEI")
  Mosquito_Female$set("public","humanEncounter",mbites_humanEncounter_manip,overwrite=TRUE)
  Mosquito_Female$set("public","timing",mbites_timing_manip,overwrite=TRUE)
  MBITES:::Parameters$set_parameters(aqua_model="emerge",
    Bs_surv=0.98,Os_surv=0.98,B_surv=0.98,O_surv=0.98,
    Bs_succeed=0.99,Os_succeed=0.99,B_succeed=0.95,O_succeed=0.99,
    S_u=0,disperse=0.1, surviveH=0.99,surviveprobeH=0.99, probeH=BASE_PROBEH,feedH=BASE_FEEDH,
    sns_a=SNS_A, sns_b=SNS_B)
}
