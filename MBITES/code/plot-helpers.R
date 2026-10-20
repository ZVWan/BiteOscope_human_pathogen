# Plotting helpers for the plot and summary scripts. The CSVs they read are written
# by the runners, one row per run (see read_complete in _common.R).

suppressMessages(library(ggplot2))

# pooled fold per sweep level, baseline = the lowest multiplier of the sweep
pooled_grid <- function(df) {
  base <- df[df$mult == min(df$mult), c("seed", "gj_vc")]
  names(base)[2] <- "b"
  m <- merge(df, base, by = "seed")
  m <- m[m$b > 0, ]
  do.call(rbind, lapply(sort(unique(m$mult)), function(mm) {
    z <- m[m$mult == mm, ]
    ci <- pooled_fold(z$gj_vc, z$b)
    data.frame(x = mm, est = ci[["est"]], lo = ci[["lo"]], hi = ci[["hi"]])
  }))
}

# the baseline group is called baseline04 in the manuscript points file and
# baseline in the DENV one
baseline_name <- function(points) {
  nm <- grep("^baseline", unique(points$scenario), value = TRUE)
  if (!length(nm)) stop("no baseline group in the points file")
  nm[1]
}

# pooled fold for one exact scenario against the baseline group, same seeds
red_point <- function(points, scen) {
  base_nm <- baseline_name(points)
  a <- points[points$scenario == scen, c("seed", "gj_vc")]
  b <- points[points$scenario == base_nm, c("seed", "gj_vc")]
  m <- merge(a, b, by = "seed")
  m <- m[m$gj_vc.y > 0, ]
  ci <- pooled_fold(m$gj_vc.x, m$gj_vc.y)
  data.frame(est = ci[["est"]], lo = ci[["lo"]], hi = ci[["hi"]], n = nrow(m))
}

# ---- devices and style -------------------------------------------------------

save_pdf <- function(obj, file, width, height) {
  force(obj)
  ok <- tryCatch({ grDevices::cairo_pdf(tempfile(), width = width, height = height)
                   grDevices::dev.off(); TRUE },
                 error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
  if (ok) grDevices::cairo_pdf(file, width = width, height = height)
  else grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off())
  print(obj)
  invisible(file)
}

ink <- "#222222"; muted <- "#6b6b6b"; scen_red <- "#CC0000"

theme_panel <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#ececec", linewidth = 0.35),
          plot.title = element_text(face = "bold", size = base_size * 1.05, colour = ink),
          plot.subtitle = element_text(size = base_size * 0.75, colour = muted, margin = margin(b = 8)),
          axis.title = element_text(size = base_size * 0.88, colour = muted),
          axis.text = element_text(colour = ink),
          plot.title.position = "plot",
          plot.margin = margin(12, 14, 10, 10))
}

# y label shared by the fold panels
fold_axis <- function() {
  scale_y_continuous("Garrett-Jones VC vs. baseline (pooled fold)",
                     labels = function(z) paste0(z, "\u00d7"))
}

# one panel: pooled fold curve over a behavioural multiplier, plus one exact point
panel_plot <- function(pg, rp, col, title, sub, xlab, label = TRUE) {
  g <- ggplot(pg, aes(x)) +
    geom_hline(yintercept = 1, linewidth = 0.4, colour = "#999999", linetype = "22") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = col, alpha = 0.14) +
    geom_line(aes(y = est), colour = col, linewidth = 1.0) +
    geom_point(aes(y = est), size = 2.4, colour = "black")
  if (!is.null(rp)) {
    g <- g +
      geom_errorbar(data = rp, aes(x = x, ymin = lo, ymax = hi), width = 0.03,
                    colour = scen_red, alpha = 0.85) +
      geom_point(data = rp, aes(x = x, y = est), size = 2.8, colour = scen_red)
    if (label)
      g <- g + geom_text(data = rp, aes(x = x, y = est, label = lab), vjust = -1.2,
                         size = 3.2, colour = scen_red)
  }
  g + scale_x_continuous(xlab, breaks = pg$x, labels = paste0(pg$x, "\u00d7")) +
    fold_axis() +
    labs(title = title, subtitle = sub) + theme_panel()
}
