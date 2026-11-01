# Infection-dependent blood-feeding behaviour.
#
# MBITES stores probeH, feedH and related parameters globally. The methods below
# are rebound to use the mosquito's pathogen state. Options hold the multipliers
# and telemetry.

library(MBITES)

# Configuration and telemetry.
# Multipliers for infectious mosquitoes; 1.0 is the baseline.
options(mbites.manip.probe_multiplier = 1.0)
options(mbites.manip.feed_multiplier  = 1.0)
# A value below 1 makes probing more costly for infectious mosquitoes.
options(mbites.manip.surviveprobe_multiplier = 1.0)
# Blood-feeding bout rate; values above 1 shorten the wait between bouts.
options(mbites.manip.bloodbout_rate_multiplier = 1.0)
# Separate multipliers for attempt (rate_b) and search (rate_bs) bouts.
options(mbites.manip.bloodattempt_rate_multiplier = 1.0)  # rate_b  (search == FALSE)
options(mbites.manip.bloodsearch_rate_multiplier  = 1.0)  # rate_bs (search == TRUE)
# Probability that an infectious mosquito survives a blood meal. This adds a
# per-meal cost because native MBITES blood meals are non-lethal.
# run_gj()/run_sim() reset it to 1.0 before each run.
set_feed_survival <- function(fs) {
  stopifnot(is.numeric(fs), length(fs) == 1L, !is.na(fs), fs > 0, fs <= 1)
  options(mbites.manip.feed_survival = fs)
  invisible(fs)
}
get_feed_survival <- function() {
  getOption("mbites.manip.feed_survival", 1.0)
}
options(mbites.manip.feed_survival = 1.0)   # default: no cost
# When TRUE, blood-feeding bouts skip flight and hazard mortality.
options(mbites.manip.risk_free_bloodbout = FALSE)
# Counters live in a mutable environment.
options(mbites.manip.telemetry = new.env(parent = emptyenv()))

reset_manip_telemetry <- function() {
  te <- getOption("mbites.manip.telemetry")
  te$probes_total      <- 0L   # human-host probes
  te$probes_infectious <- 0L   # probes by infectious mosquitoes
  te$feeds_infectious  <- 0L   # full meals by infectious mosquitoes
  te$probedeaths_inf   <- 0L   # infectious probe deaths
  te$feeddeaths_inf    <- 0L   # infectious feed-cost deaths
  te$bloodbouts_inf    <- 0L   # blood bouts by infectious mosquitoes
  te$attemptbouts_inf  <- 0L   # blood attempt bouts by infectious mosquitoes
  invisible(te)
}
reset_manip_telemetry()

# Expose the pathogen's infectious state.
SEI_Pathogen$set(
  which = "public", name = "is_infectious",
  value = function() private$infectious, overwrite = TRUE
)

# Return TRUE when the mosquito has an infectious pathogen.
# The incubating stage does not count.
mbites_is_infectious_now <- function() {
  for (p in private$pathogens) {
    if ("SEI_Pathogen" %in% class(p) && isTRUE(p$is_infectious())) return(TRUE)
  }
  FALSE
}
Mosquito_Female$set(
  which = "public", name = "is_infectious_now",
  value = mbites_is_infectious_now, overwrite = TRUE
)

# Human encounter, with infection-dependent multipliers and counters.
mbites_humanEncounter_manip <- function() {
  rf <- getOption("mbites.manip.risk_free_bloodbout", FALSE)
  if (!rf && runif(1) < 1 - MBITES:::Parameters$get_surviveH()) {
    # dies before probing
    private$alive <- FALSE
    private$cod   <- "surviveH"
  } else {
    # Apply infection-dependent multipliers.
    probeH        <- MBITES:::Parameters$get_probeH()
    feedH         <- MBITES:::Parameters$get_feedH()
    surviveprobeH <- MBITES:::Parameters$get_surviveprobeH()
    infnow <- self$is_infectious_now()
    if (infnow) {
      probeH        <- min(1, probeH * getOption("mbites.manip.probe_multiplier", 1.0))
      feedH         <- min(1, feedH  * getOption("mbites.manip.feed_multiplier",  1.0))
      surviveprobeH <- min(1, surviveprobeH *
                              getOption("mbites.manip.surviveprobe_multiplier", 1.0))
    }
    if (runif(1) < probeH) {
      # Count the probe.
      te <- getOption("mbites.manip.telemetry")
      te$probes_total <- te$probes_total + 1L
      if (infnow) te$probes_infectious <- te$probes_infectious + 1L

      # Probe the host.
      self$trackProbe()
      self$probeHost()   # PATHOGEN-XX.R (mosquito -> human transmission)
      if (!rf && runif(1) < 1 - surviveprobeH) {
        # dies before feeding
        if (infnow) te$probedeaths_inf <- te$probedeaths_inf + 1L
        private$alive <- FALSE
        private$cod   <- "surviveprobeH"
      } else {
        if (runif(1) < feedH) {
          # Start blood feeding.
          if (infnow) te$feeds_infectious <- te$feeds_infectious + 1L
          self$trackFeed()
          self$feedHost()   # human -> mosquito transmission
          self$BloodMeal()
          # Extra cost during engorgement.
          if (infnow) {
            fs <- get_feed_survival()
            if (fs < 1 && runif(1) < 1 - fs) {
              private$alive <- FALSE
              private$cod   <- "manip_feedcost"
              te$feeddeaths_inf <- te$feeddeaths_inf + 1L
            }
          }
        }
      }
    }
  }
}
Mosquito_Female$set(
  which = "public", name = "humanEncounter",
  value = mbites_humanEncounter_manip, overwrite = TRUE
)

# Timing; infectious mosquitoes can start the next blood bout sooner.
mbites_timing_manip <- function() {
  if (private$search) {
    private$tNext <- switch(private$state,
      B = MBITES:::Parameters$ttEvent$BoutBs(private$tNow),
      O = MBITES:::Parameters$ttEvent$BoutOs(private$tNow),
      M = MBITES:::Parameters$ttEvent$BoutMs(private$tNow),
      S = MBITES:::Parameters$ttEvent$BoutSs(private$tNow))
  } else {
    private$tNext <- switch(private$state,
      B = MBITES:::Parameters$ttEvent$BoutB(private$tNow),
      O = MBITES:::Parameters$ttEvent$BoutO(private$tNow),
      M = MBITES:::Parameters$ttEvent$BoutM(private$tNow),
      S = MBITES:::Parameters$ttEvent$BoutS(private$tNow))
  }
  # Scale the waiting time for infectious mosquitoes.
  if (identical(private$state, "B") && self$is_infectious_now()) {
    k_all <- getOption("mbites.manip.bloodbout_rate_multiplier", 1.0)
    k_sub <- if (isTRUE(private$search)) getOption("mbites.manip.bloodsearch_rate_multiplier", 1.0)
             else                        getOption("mbites.manip.bloodattempt_rate_multiplier", 1.0)
    k <- k_all * k_sub
    if (k != 1.0) private$tNext <- private$tNow + (private$tNext - private$tNow) / k
    te <- getOption("mbites.manip.telemetry")
    te$bloodbouts_inf <- te$bloodbouts_inf + 1L
    if (!isTRUE(private$search)) te$attemptbouts_inf <- te$attemptbouts_inf + 1L
  }
}
Mosquito_Female$set(
  which = "public", name = "timing",
  value = mbites_timing_manip, overwrite = TRUE
)

# Skip flight and hazard mortality during blood-feeding bouts when requested.
mbites_survival_manip <- function() {
  if (getOption("mbites.manip.risk_free_bloodbout", FALSE) && identical(private$state, "B")) {
    return(invisible())   # risk-free: skip mortality on blood bouts
  }
  self$surviveFlight()
  self$surviveHazards()
}
Mosquito_Female$set(which = "public", name = "survival",
                    value = mbites_survival_manip, overwrite = TRUE)

# Keep the timing update and trackRest, but skip the death draw when risk-free.
mbites_PPRFlight_manip <- function() {
  if (private$alive) {
    private$tNow <- MBITES:::Parameters$ttEvent$ppr(private$tNow)
    self$trackRest()
    if (!getOption("mbites.manip.risk_free_bloodbout", FALSE)) {
      if (runif(1) < 1 - self$pPPRFlight()) {
        private$alive <- FALSE
        private$cod   <- "PPRFlight"
      }
    }
  }
}
Mosquito_Female$set(which = "public", name = "PPRFlight",
                    value = mbites_PPRFlight_manip, overwrite = TRUE)

cat("[manip] infection-conditional probe/feed/timing/risk-free manipulation installed on Mosquito_Female\n")
