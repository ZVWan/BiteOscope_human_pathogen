# The probing scenario under both pathogen settings in one figure:
#   malaria settings  results/gj-probing-30seeds.csv     (30seeds_sweeps/)
#   dengue settings   results/gj-denv-probing-30seeds.csv (denv_sweeps/)
#   red dots          S3_probe and S4_probe (malaria settings)
#                     S4_probe (dengue settings)
#
# Each curve is the pooled fold against its own no-manipulation baseline, so both sit
# on the same axis: 1 means "unchanged against that arm's own baseline". The two arms
# are separate simulations (same seeds, but the parameter change sends the RNG streams
# apart), so read the difference descriptively, not as a paired test.
#
# Usage: Rscript plot_probing_malaria_dengue.R   (needs all four runners done)
# Output: 30seeds_sweeps/results/probing-malaria-dengue.pdf

COLOUR_MALARIA <- "#009E73"    # green, the colour of the manuscript probing panel
COLOUR_DENGUE  <- "#E69F00"    # amber: pure yellow loses the line against white
COLOUR_POINT   <- "#CC0000"    # both exact points stay red

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
OUT_DENV <- file.path(ROOT, "denv_sweeps", "results")
suppressMessages(source(file.path(ROOT, "code", "_common.R")))
suppressMessages(source(file.path(ROOT, "code", "plot-helpers.R")))

mal <- pooled_grid(read_complete(file.path(OUT, "gj-probing-30seeds.csv")))
den <- pooled_grid(read_complete(file.path(OUT_DENV, "gj-denv-probing-30seeds.csv")))
if (!identical(mal$x, den$x))
  stop("the two sweeps were run on different multiplier grids: ",
       paste(mal$x, collapse = " "), " vs ", paste(den$x, collapse = " "))
if (!isTRUE(all.equal(mal$est[1], 1)) || !isTRUE(all.equal(den$est[1], 1)))
  stop("a curve does not pass through 1 at x = 1, so the baselines are not the lowest multiplier")

pmal <- read_complete(file.path(OUT, "gj-exact-points-30seeds.csv"))
pden <- read_complete(file.path(OUT_DENV, "gj-denv-points-30seeds.csv"))
rp <- rbind(
  cbind(red_point(pmal, "S3_probe"), x = 1.29, lab = "S3  NF135, Anopheles",   ctx = "malaria"),
  cbind(red_point(pmal, "S4_probe"), x = 1.52, lab = "S4  malaria settings",   ctx = "malaria"),
  cbind(red_point(pden, "S4_probe"), x = 1.52, lab = "S4  dengue settings",    ctx = "dengue"))
# where two points share a multiplier, put the higher label above and the lower below
rp$vj <- ave(rp$est, rp$x, FUN = function(v) ifelse(rank(v) == length(v), -1.6, 2.0))

curves <- rbind(cbind(mal, ctx = "malaria settings"),
                cbind(den, ctx = "dengue settings"))
cols <- c("malaria settings" = COLOUR_MALARIA, "dengue settings" = COLOUR_DENGUE)

g <- ggplot(curves, aes(x, est, colour = ctx, fill = ctx)) +
  geom_hline(yintercept = 1, linewidth = 0.4, colour = "#999999", linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.9, colour = "black", show.legend = FALSE) +
  geom_errorbar(data = rp, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0.03, colour = COLOUR_POINT, alpha = 0.85) +
  geom_point(data = rp, aes(x = x, y = est), inherit.aes = FALSE,
             size = 2.8, colour = COLOUR_POINT) +
  geom_text(data = rp, aes(x = x, y = est, label = lab, vjust = vj), inherit.aes = FALSE,
            size = 3.0, colour = COLOUR_POINT) +
  scale_x_continuous("probability of probing (probeH, \u00d7 baseline)",
                     breaks = mal$x, labels = paste0(mal$x, "\u00d7")) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  fold_axis() +
  labs(title = "Probing and Garrett-Jones VC under malaria and dengue pathogen settings (pooled fold)",
       subtitle = paste0("EIP ", EIP_DAYS, " d in model and readout; each curve is the pooled fold against ",
                         "its own no-manipulation baseline over ", rp$n[1],
                         " seeds; band = paired bootstrap 95% CI; red = exact scenario points."),
       colour = NULL, fill = NULL) +
  theme_panel() +
  theme(legend.position = "top", legend.title = element_blank(), legend.margin = margin(b = -4))

f <- file.path(OUT, "probing-malaria-dengue.pdf")
save_pdf(g, f, width = 5.33, height = 3.8)

cat(sprintf("probeH sweep at EIP %d d (fold against each arm's own baseline):\n", EIP_DAYS))
cat(sprintf("  %-6s  %-22s  %-22s\n", "x", "malaria settings", "dengue settings"))
for (i in seq_len(nrow(mal)))
  cat(sprintf("  %-6.2f  %.3f [%.3f, %.3f]      %.3f [%.3f, %.3f]\n",
              mal$x[i], mal$est[i], mal$lo[i], mal$hi[i], den$est[i], den$lo[i], den$hi[i]))
cat("\nexact points (red):\n")
for (i in seq_len(nrow(rp)))
  cat(sprintf("  %-24s x = %.2f  %.3f [%.3f, %.3f]  n = %d\n",
              rp$lab[i], rp$x[i], rp$est[i], rp$lo[i], rp$hi[i], rp$n[i]))
cat("wrote", f, "\n")
