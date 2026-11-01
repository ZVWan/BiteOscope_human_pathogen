# Per-mosquito counters for the "infectious probes per infected mosquito" readout.
# Sourced after MBITES/code/_common.R (or denv-pathogen.R); the package and the
# model code are not changed, two Mosquito_Female methods are rebound with the
# package body plus counting:
#   feedHost  records the time a mosquito first acquires the pathogen
#   probeHost counts human probes by infectious mosquitoes, and the probes that
#             actually infected the human, per mosquito
# No random numbers are drawn here, so the simulation itself is unchanged.

reset_metric <- function() {
  options(mbites.metric = list(tinf   = new.env(hash = TRUE),
                               probes = new.env(hash = TRUE),
                               tx     = new.env(hash = TRUE)))
}
reset_metric()

.bump <- function(e, k) assign(k, get0(k, envir = e, inherits = FALSE, ifnotfound = 0L) + 1L,
                               envir = e)

metric_feedHost <- function() {
  had <- length(private$pathogens) > 0
  host <- MBITES:::Globals$get_tile(private$tileID)$get_human(private$hostID)
  to_add <- host$feedHost(self, private$pathogens)
  for (mPathogen in to_add) if (!is.null(mPathogen)) self$add_pathogen(mPathogen)
  if (!had && length(private$pathogens) > 0)
    assign(as.character(self$get_id()), private$tNow, envir = getOption("mbites.metric")$tinf)
}
Mosquito_Female$set(which = "public", name = "feedHost", value = metric_feedHost, overwrite = TRUE)

metric_probeHost <- function() {
  host <- MBITES:::Globals$get_tile(private$tileID)$get_human(private$hostID)
  if (is.null(host)) stop(paste("the host is null", private$hostID))
  is_SEI <- function(p) "SEI_Pathogen" %in% class(p)
  inf <- self$is_infectious_now()
  if (inf) {
    mt <- getOption("mbites.metric"); k <- as.character(self$get_id())
    .bump(mt$probes, k)
    h_before <- any(as.logical(host$on_pathogens(is_SEI)))
  }
  for (pathogen in private$pathogens) pathogen$probeHost(self, host)
  host$mosquitoProbe(self$get_id())
  if (inf && !h_before && any(as.logical(host$on_pathogens(is_SEI)))) .bump(mt$tx, k)
}
Mosquito_Female$set(which = "public", name = "probeHost", value = metric_probeHost, overwrite = TRUE)

# Summary per run: mosquitoes infected by day `cut`, their infectious probes and
# the human infections they caused (cut-offs leave time to live out the infection)
metric_summary <- function(cuts = c(60, 70)) {
  mt <- getOption("mbites.metric")
  ids <- ls(mt$tinf); tinf <- vapply(ids, function(k) get(k, envir = mt$tinf), numeric(1))
  cnt <- function(e) vapply(ids, function(k) get0(k, envir = e, inherits = FALSE, ifnotfound = 0L),
                            integer(1))
  pr <- cnt(mt$probes); tx <- cnt(mt$tx)
  out <- list()
  for (c in cuts) {
    sel <- tinf <= c
    out[[paste0("ninf_", c)]]   <- sum(sel)
    out[[paste0("probes_", c)]] <- sum(pr[sel])
    out[[paste0("tx_", c)]]     <- sum(tx[sel])
  }
  as.data.frame(out)
}
