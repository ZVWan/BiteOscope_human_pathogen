# The two DENV panels from the 30-seed CSVs: pooled fold against the
# behavioural multiplier, with the DENV exact point in red.
#
# Usage: Rscript plot_denv_panels.R   (run the two DENV runners first)
# Outputs: results/denv-{rateb,probing}-pooled.pdf

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
suppressMessages(source(file.path(ROOT, "code", "denv-pathogen.R")))
suppressMessages(source(file.path(ROOT, "code", "plot-helpers.R")))

points <- read_complete(file.path(OUT, "gj-denv-points-30seeds.csv"))
sub_common <- sprintf("DENV pathogen settings, EIP %d d in model and readout; pooled fold = %s; band = paired bootstrap 95%% CI.",
                      EIP_DAYS, "\u03a3 VC(\u00d7m)/\u03a3 VC(\u00d71) over 30 seeds")

panels <- list(
  list(csv = "gj-denv-rateb-30seeds.csv", stem = "denv-rateb-pooled", col = "#D55E00",
       scen = "S2_rate", x = 1.6, lab = "S2 (DENV)",
       title = "DENV landing rate and Garrett-Jones VC (pooled fold)",
       xlab = "landing rate of infectious mosquitoes (rate_b, \u00d7 baseline)"),
  list(csv = "gj-denv-probing-30seeds.csv", stem = "denv-probing-pooled", col = "#0072B2",
       scen = "S4_probe", x = 1.52, lab = "S4 (DENV)",
       title = "DENV probing and Garrett-Jones VC (pooled fold)",
       xlab = "probability of probing (probeH, \u00d7 baseline)")
)

for (p in panels) {
  pg <- pooled_grid(read_complete(file.path(OUT, p$csv)))
  rp <- red_point(points, p$scen); rp$x <- p$x; rp$lab <- p$lab
  sub <- sprintf("%s Red: pooled %s (n = %d).", sub_common, p$lab, rp$n)
  g <- panel_plot(pg, rp, p$col, p$title, sub, p$xlab)
  f <- file.path(OUT, paste0(p$stem, ".pdf"))
  save_pdf(g, f, width = 5.33, height = 3.8)
  cat(" ", p$stem, "\n")
  for (i in seq_len(nrow(pg)))
    cat(sprintf("    sweep x=%.2f  %.3f [%.3f, %.3f]\n", pg$x[i], pg$est[i], pg$lo[i], pg$hi[i]))
  cat(sprintf("    %-12s %.3f [%.3f, %.3f]\n", p$lab, rp$est, rp$lo, rp$hi))
  cat("    wrote", f, "\n")
}
