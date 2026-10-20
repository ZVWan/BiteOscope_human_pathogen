# Fig. 4c variant: the engorgement sweep drawn on the same x axis as the other
# panels (1x to 2.5x) but with the pooled curve stopping at 2x. feedH is 0.5 at
# baseline, so it reaches its ceiling of 1 at 2x: the 2.25x and 2.5x multipliers
# repeat the 2x configuration and add nothing but a flat tail.
#
# Usage: Rscript plot_feed_panel_2x.R   (run run_sweeps.R and run_exact_points.R first)
# Outputs: results/feed-panel-2x.pdf

.fp_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (length(p) && nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
SCRIPT_DIR <- .fp_dir()
ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))
OUT <- file.path(SCRIPT_DIR, "results")
suppressMessages(source(file.path(ROOT, "code", "_common.R")))
suppressMessages(source(file.path(ROOT, "code", "plot-helpers.R")))

X_MAX_DATA <- 2      # the pooled curve stops here
X_MAX_AXIS <- 2.5    # the axis keeps the span of the other panels
ticks <- seq(1, X_MAX_AXIS, by = 0.25)

points <- read_complete(file.path(OUT, "gj-exact-points-30seeds.csv"))
pg <- pooled_grid(read_complete(file.path(OUT, "gj-feed-30seeds.csv")))
pg <- pg[pg$x <= X_MAX_DATA, ]
rp <- red_point(points, "S5_feed"); rp$x <- 1.92; rp$lab <- "S5 (NF54)"

sub_common <- paste0("EIP ", EIP_DAYS, " d in model and readout; pooled fold = \u03a3 VC(\u00d7m)/",
                     "\u03a3 VC(\u00d71) over 30 seeds; band = paired bootstrap 95% CI.")
sub <- sprintf("%s Red: pooled %s (n = %s).", sub_common, rp$lab, rp$n)
xlab <- "probability of engorgement (feedH, \u00d7 baseline)"

g <- suppressMessages(
  panel_plot(pg, rp, "#009E73",
             "Probability of engorgement and Garrett-Jones VC (pooled fold)",
             sub, xlab) +
  # panel_plot() sets the x scale from the plotted multipliers; this one keeps the
  # axis span of the other panels while the curve itself ends at 2x
  scale_x_continuous(xlab, limits = c(1, X_MAX_AXIS), breaks = ticks,
                     labels = paste0(ticks, "\u00d7")))

f <- file.path(OUT, "feed-panel-2x.pdf")
suppressMessages(save_pdf(g, f, width = 5.33, height = 3.8))

cat(sprintf("panel values at EIP %d d (curve to x%.2f, axis to x%.2f):\n",
            EIP_DAYS, X_MAX_DATA, X_MAX_AXIS))
for (i in seq_len(nrow(pg)))
  cat(sprintf("    sweep x=%.2f  %.3f [%.3f, %.3f]\n", pg$x[i], pg$est[i], pg$lo[i], pg$hi[i]))
cat(sprintf("    %-12s %.3f [%.3f, %.3f]\n", rp$lab, rp$est, rp$lo, rp$hi))
cat("    wrote", f, "\n")
