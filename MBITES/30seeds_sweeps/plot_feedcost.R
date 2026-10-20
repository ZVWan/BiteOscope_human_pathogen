# Panel d from the S6 CSV: the engorgement boost with and without a
# per-meal survival cost, both against a no-manipulation baseline.
#
# Usage: Rscript plot_feedcost.R     (run run_feedcost.R first)
# Output: results/feedcost-pooled.pdf

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

d <- read_complete(file.path(OUT, "gj-feedcost-30seeds.csv"))

fold_at_fs <- function(scen_name, anchor_rows) {
  base <- anchor_rows[, c("seed", "gj_vc")]; names(base)[2] <- "b"
  do.call(rbind, lapply(sort(unique(d$feed_survival[d$scenario == scen_name])), function(fs) {
    a <- d[d$scenario == scen_name & abs(d$feed_survival - fs) < 1e-9, c("seed", "gj_vc")]
    m <- merge(a, base, by = "seed"); m <- m[m$b > 0, ]
    ci <- pooled_fold(m$gj_vc, m$b)
    data.frame(death = 1 - fs, fs = fs, est = ci[["est"]], lo = ci[["lo"]], hi = ci[["hi"]])
  }))
}

rb <- fold_at_fs("feedcost", d[d$scenario == "baseline", ])
rn <- fold_at_fs("feedcost_noboost", d[d$scenario == "feedcost_noboost" & abs(d$feed_survival - 1) < 1e-9, ])
rb$arm <- "feedH \u00d71.92 (engorgement boost)"
rn$arm <- "feedH \u00d71 (behaviour unchanged)"
p <- rbind(rb, rn)
p$arm <- factor(p$arm, levels = c("feedH \u00d71.92 (engorgement boost)",
                                  "feedH \u00d71 (behaviour unchanged)"))

anchor <- rn$est[abs(rn$fs - 1) < 1e-9]
stopifnot(length(anchor) == 1, abs(anchor - 1) < 1e-9)   # no-boost arm at fs = 1.00 is the anchor

cols <- c("feedH \u00d71.92 (engorgement boost)" = "#009E73",
          "feedH \u00d71 (behaviour unchanged)" = "#0072B2")
g <- ggplot(p, aes(death, est, colour = arm, fill = arm)) +
  geom_hline(yintercept = 1, linewidth = 0.4, colour = "#999999", linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2, colour = "black") +
  annotate("text", x = max(p$death) * 1.02, y = 1.03, label = "no net gain (fold = 1)",
           colour = muted, size = 3, hjust = 1, vjust = 0) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  scale_x_continuous("per-meal host-defensive mortality of infectious mosquitoes (1 \u2212 feed_survival)",
                     breaks = seq(0, max(p$death), by = 0.05)) +
  fold_axis() +
  labs(title = "S5 engorgement scenario under per-meal survival cost, with and without the behavioural boost",
       subtitle = paste0("EIP ", EIP_DAYS, " d in model and readout; pooled fold vs no-manipulation baseline; ",
                         "bands = paired bootstrap 95% CI; the gap between the arms at a given cost is the ",
                         "contribution of feedH \u00d71.92."),
       colour = NULL, fill = NULL) +
  theme_panel() +
  theme(legend.position = "top", legend.title = element_blank(), legend.margin = margin(b = -4))

f <- file.path(OUT, "feedcost-pooled.pdf")
save_pdf(g, f, width = 5.33, height = 3.8)

cat(sprintf("panel d values at EIP %d d (fs = 1 - death):\n", EIP_DAYS))
for (i in seq_len(nrow(p)))
  cat(sprintf("  %-34s death %.2f  fs %.2f  %.3f [%.3f, %.3f]\n",
              as.character(p$arm[i]), p$death[i], p$fs[i], p$est[i], p$lo[i], p$hi[i]))
cat("wrote", f, "\n")
