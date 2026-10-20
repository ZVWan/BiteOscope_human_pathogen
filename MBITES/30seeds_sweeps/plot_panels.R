# The three sweep panels from the 30-seed CSVs: pooled fold against
# the behavioural multiplier, with the exact scenario points in red.
#
#   a  landing rate (rate_b)   red S2 x1.6    DENV
#   b  probing (probeH)        red S3 x1.29 and S4 x1.52   NF135, DENV
#   c  engorgement (feedH)     red S5 x1.92   NF54, drawn to x2 only
#
# Panel c stops at x2 because feedH is 0.5 at baseline, so the x2.25 and x2.5
# multipliers both land on the ceiling of 1 and repeat the x2 run exactly.
#
# Usage: Rscript plot_panels.R     (run run_sweeps.R and run_exact_points.R first)
# Outputs: results/{rateb,probing,feed}-pooled.pdf

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
suppressMessages(source(file.path(ROOT, "code", "_common.R")))
suppressMessages(source(file.path(ROOT, "code", "plot-helpers.R")))

points <- read_complete(file.path(OUT, "gj-exact-points-30seeds.csv"))
sub_common <- paste0("EIP ", EIP_DAYS, " d in model and readout; pooled fold = ",
                     "\u03a3 VC(\u00d7m)/\u03a3 VC(\u00d71) over 30 seeds; band = paired bootstrap 95% CI.")

panels <- list(
  list(csv = "gj-rateb-30seeds.csv", stem = "rateb-pooled", col = "#D55E00",
       scen = "S2_rate04", x = 1.6, lab = "S2 (DENV)", xmax = 2.5,
       title = "Higher landing rate LOWERS Garrett-Jones VC at the common probing baseline (pooled fold)",
       xlab = "landing rate of infectious mosquitoes (rate_b, \u00d7 baseline)"),
  list(csv = "gj-probing-30seeds.csv", stem = "probing-pooled", col = "#0072B2",
       scen = c("S3_probe", "S4_probe"), x = c(1.29, 1.52), lab = c("S3 (NF135)", "S4 (DENV)"),
       xmax = 2.5,
       title = "Enhanced probing raises Garrett-Jones VC (pooled fold)",
       xlab = "probability of probing (probeH, \u00d7 baseline)"),
  list(csv = "gj-feed-30seeds.csv", stem = "feed-pooled", col = "#009E73",
       scen = "S5_feed", x = 1.92, lab = "S5 (NF54)", xmax = 2,
       title = "Probability of engorgement and Garrett-Jones VC (pooled fold)",
       xlab = "probability of engorgement (feedH, \u00d7 baseline)")
)

cat(sprintf("panel values at EIP %d d:\n", EIP_DAYS))
for (p in panels) {
  pg <- pooled_grid(read_complete(file.path(OUT, p$csv)))
  pg <- pg[pg$x <= p$xmax, ]
  rp <- do.call(rbind, lapply(seq_along(p$scen), function(i) {
    r <- red_point(points, p$scen[i]); r$x <- p$x[i]; r$lab <- p$lab[i]; r
  }))
  sub <- sprintf("%s Red: pooled %s (n = %s).", sub_common,
                 paste(p$lab, collapse = ", "), paste(rp$n, collapse = "/"))
  g <- panel_plot(pg, rp, p$col, p$title, sub, p$xlab)
  f <- file.path(OUT, paste0(p$stem, ".pdf"))
  save_pdf(g, f, width = 5.33, height = 3.8)
  cat(" ", p$stem, "\n")
  for (i in seq_len(nrow(pg)))
    cat(sprintf("    sweep x=%.2f  %.3f [%.3f, %.3f]\n", pg$x[i], pg$est[i], pg$lo[i], pg$hi[i]))
  for (i in seq_len(nrow(rp)))
    cat(sprintf("    %-12s %.3f [%.3f, %.3f]\n", rp$lab[i], rp$est[i], rp$lo[i], rp$hi[i]))
  cat("    wrote", f, "\n")
}
