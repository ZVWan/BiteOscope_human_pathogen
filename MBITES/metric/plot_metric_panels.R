# Fig. 4 panels with the per-infected-female readout.
#
#   readout   infectious probes per infected female (females infected by day CUT),
#             pooled over seeds: sum(probes) / sum(infected females)
#   fold      manipulated / baseline on the same seeds; band = paired percentile
#             bootstrap 95% CI over seeds (B = 2000)
#
# Usage: Rscript plot_metric_panels.R   (works on partial results, too)
# Outputs in results/: fig4{a,b,c,d}-metric.pdf, fig4-metric.pdf, metric-values.csv

suppressMessages(library(ggplot2)); library(grid)
options(warn = 1)

# Paths are resolved from this script, so the package runs from any checkout.
# MBITES_MODEL_DIR can point at another copy of the model tree.
SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (length(p) && nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
})
HERE <- normalizePath(SCRIPT_DIR)
EN   <- normalizePath(Sys.getenv("MBITES_MODEL_DIR", file.path(HERE, "..", "MBITES")))
OUT  <- Sys.getenv("METRIC_OUT", file.path(HERE, "results"))
source(file.path(HERE, "plot-helpers.R"))
CUT <- 70
XMAX <- 2      # sweeps are plotted up to this multiplier (runs go to 2.5)
NUM <- paste0("probes_", CUT); DEN <- paste0("ninf_", CUT)

sw <- read.csv(file.path(OUT, "metric-sweeps-30seeds.csv"), stringsAsFactors = FALSE)
pt <- read.csv(file.path(OUT, "metric-points-30seeds.csv"), stringsAsFactors = FALSE)

# pooled fold of manipulated rows m against baseline rows b, paired by seed
fold <- function(m, b, B = 2000) {
  s <- intersect(m$seed, b$seed)
  m <- m[match(s, m$seed), ]; b <- b[match(s, b$seed), ]
  r <- function(i) (sum(m[i, NUM]) / sum(m[i, DEN])) / (sum(b[i, NUM]) / sum(b[i, DEN]))
  set.seed(20261030)
  bs <- replicate(B, r(sample(length(s), replace = TRUE)))
  data.frame(est = r(seq_along(s)), lo = unname(quantile(bs, .025)),
             hi = unname(quantile(bs, .975)), n = length(s))
}

sweep_curve <- function(set, study) {
  z <- sw[sw$set == set & sw$study == study & sw$x <= XMAX + 1e-9, ]
  if (!nrow(z)) return(NULL)
  b <- z[z$x == 1, ]
  do.call(rbind, lapply(sort(unique(z$x)), function(x)
    cbind(x = x, fold(z[z$x == x, ], b))))
}
point <- function(set, scen, x, lab) {
  z <- pt[pt$set == set, ]
  cbind(fold(z[z$scenario == scen, ], z[grepl("^baseline", z$scenario), ]),
        x = x, lab = lab)
}

y_new <- scale_y_continuous("infectious probes per\ninfected female (fold)",
                            labels = function(z) paste0(z, "\u00d7"))
sub_new <- sprintf("Infected by day %d; 30 seeds; band = bootstrap 95%% CI; red = exact scenarios", CUT)
COL <- c(malaria = "#0072B2", dengue = "#E69F00")

# a, c: one curve and its exact point
single_panel <- function(set, study, scen, x, lab, col, title, xlab) {
  pg <- sweep_curve(set, study); if (is.null(pg)) return(NULL)
  rp <- point(set, scen, x, lab)
  suppressMessages(panel_plot(pg, rp, col, title, sub_new, xlab) + y_new)
}
pa <- single_panel("dengue", "rateb", "S2_rate", 1.6, "S2 (DENV)", "#D55E00",
                   "a  Attempt rate (dengue settings)",
                   "blood-feeding attempt rate (rate_b, \u00d7 baseline)")
pc <- single_panel("malaria", "feed", "S5_feed", 1.92, "S5 (NF54)", "#009E73",
                   "c  Engorgement (malaria settings)",
                   "probability of engorgement (feedH, \u00d7 baseline)")

# b: probing under both pathogen settings
lab_curve <- function(r, ctx) if (is.null(r)) NULL else cbind(r, ctx = ctx)
curves_b <- rbind(lab_curve(sweep_curve("malaria", "probing"), "malaria settings"),
                  lab_curve(sweep_curve("dengue", "probing"), "dengue settings"))
rp_b <- rbind(point("malaria", "S3_probe", 1.29, "S3  NF135"),
              point("malaria", "S4_probe", 1.52, "S4  malaria settings"),
              point("dengue", "S4_probe", 1.52, "S4  dengue settings"))
# labels: S3 to the left, S4 malaria to the lower right, S4 dengue to the upper left
rp_b$vj <- c(0.5, 1.4, -0.6); rp_b$hj <- c(1.12, -0.12, 1.12)
cols_b <- c("malaria settings" = unname(COL["malaria"]), "dengue settings" = unname(COL["dengue"]))
pb <- if (!is.null(curves_b)) ggplot(curves_b, aes(x, est, colour = ctx, fill = ctx)) +
  geom_hline(yintercept = 1, linewidth = 0.4, colour = "#999999", linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.0) + geom_point(size = 1.9, colour = "black", show.legend = FALSE) +
  geom_errorbar(data = rp_b, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0.03, colour = scen_red, alpha = 0.85) +
  geom_point(data = rp_b, aes(x = x, y = est), inherit.aes = FALSE, size = 2.8, colour = scen_red) +
  geom_text(data = rp_b, aes(x = x, y = est, label = lab, vjust = vj, hjust = hj), inherit.aes = FALSE,
            size = 3.0, colour = scen_red) +
  scale_x_continuous("probability of probing (probeH, \u00d7 baseline)",
                     breaks = sort(unique(curves_b$x)), labels = paste0(sort(unique(curves_b$x)), "\u00d7")) +
  scale_colour_manual(values = cols_b) + scale_fill_manual(values = cols_b) + y_new +
  labs(title = "b  Probing (malaria and dengue settings)", subtitle = sub_new,
       colour = NULL, fill = NULL) +
  theme_panel() + theme(legend.position = "top", legend.margin = margin(b = -4))

# d: S6 feed cost, boost arm vs no-manipulation baseline, no-boost arm vs its own fs = 1
fcd <- sw[sw$study == "feedcost", ]
arm_curve <- function(scen, base)
  do.call(rbind, lapply(sort(unique(fcd$x[fcd$scenario == scen])), function(fs)
    cbind(death = 1 - fs, fold(fcd[fcd$scenario == scen & abs(fcd$x - fs) < 1e-9, ], base))))
d_arms <- function() {
  if (!any(fcd$scenario == "baseline")) return(NULL)
  rb <- arm_curve("feedcost", fcd[fcd$scenario == "baseline", ])
  rn <- arm_curve("feedcost_noboost",
                  fcd[fcd$scenario == "feedcost_noboost" & abs(fcd$x - 1) < 1e-9, ])
  rbind(cbind(rb, arm = "feedH \u00d71.92 (engorgement boost)"),
        cbind(rn, arm = "feedH \u00d71 (behaviour unchanged)"))
}
cols_d <- c("feedH \u00d71.92 (engorgement boost)" = "#009E73",
            "feedH \u00d71 (behaviour unchanged)" = "#0072B2")
d_plot <- function(p, y, title, sub)
  ggplot(p, aes(death, est, colour = arm, fill = arm)) +
  geom_hline(yintercept = 1, linewidth = 0.4, colour = "#999999", linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2, colour = "black") +
  scale_colour_manual(values = cols_d) + scale_fill_manual(values = cols_d) +
  scale_x_continuous("per-meal mortality (1 - feed_survival)",
                     breaks = seq(0, 0.3, by = 0.05)) + y +
  labs(title = title, subtitle = sub, colour = NULL, fill = NULL) +
  theme_panel() + theme(legend.position = "top", legend.margin = margin(b = -4)) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
p_d <- d_arms()
pd <- if (!is.null(p_d)) d_plot(p_d, y_new, "d  Engorgement boost + per-meal cost (malaria)",
                                sub_new)

panels <- list(a = pa, b = pb, c = pc, d = pd)
for (k in names(panels)) if (!is.null(panels[[k]]))
  save_pdf(panels[[k]], file.path(OUT, sprintf("fig4%s-metric.pdf", k)), width = 5.33, height = 3.8)

grid_page <- function(plots, file, ncol = 2, w = 5.33, h = 3.8) {
  plots <- Filter(Negate(is.null), plots); nrow <- ceiling(length(plots) / ncol)
  grDevices::pdf(file, width = w * ncol, height = h * nrow)
  grid.newpage(); pushViewport(viewport(layout = grid.layout(nrow, ncol)))
  for (i in seq_along(plots))
    print(plots[[i]], vp = viewport(layout.pos.row = (i - 1) %/% ncol + 1,
                                    layout.pos.col = (i - 1) %% ncol + 1))
  grDevices::dev.off(); invisible(file)
}
grid_page(panels, file.path(OUT, "fig4-metric.pdf"))

# ---- values table ----------------------------------------------------------------
vals <- list()
for (s in list(c("dengue", "rateb"), c("malaria", "rateb"), c("malaria", "probing"),
               c("dengue", "probing"), c("malaria", "feed"))) {
  r <- sweep_curve(s[1], s[2]); if (!is.null(r))
    vals[[length(vals) + 1]] <- cbind(series = paste(s[2], s[1]), r)
}
for (p in list(c("dengue", "S2_rate"), c("malaria", "S2_rate04"), c("malaria", "S3_probe"),
               c("malaria", "S4_probe"), c("dengue", "S4_probe"), c("malaria", "S5_feed")))
  vals[[length(vals) + 1]] <- cbind(series = paste(p[2], p[1]), point(p[1], p[2], NA, "")[, 1:5])
if (!is.null(p_d)) { r <- d_arms()
  vals[[length(vals) + 1]] <- data.frame(series = paste("feedcost", r$arm), x = r$death,
                                         est = r$est, lo = r$lo, hi = r$hi, n = r$n) }
vals <- do.call(rbind, lapply(vals, function(v) v[, c("series", "x", "est", "lo", "hi", "n")]))
write.csv(vals, file.path(OUT, "metric-values.csv"), row.names = FALSE)
cat("seeds per sweep level (min):", min(vals$n), "\n")
cat("wrote figures and metric-values.csv to", OUT, "\n")
