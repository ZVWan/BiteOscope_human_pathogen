# Every value the Fig. 4 panels draw, collected into one table: sweep curves, the
# exact scenario points and the S6 arms. Each scenario is taken from its own
# parameter set - the dengue virus runs for the landing-rate scenario and for the
# DENV1 probing point, the Plasmodium runs for the rest - so the table matches the
# figure. Reads the CSVs in the results/ folders only.
#
# Usage: Rscript summary_values.R
# Output: results/gj-values-30seeds.csv and the console tables.

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
SEEDS_EXPECTED <- 30L        # the seed count the runners use
suppressMessages(source(file.path(ROOT, "code", "_common.R")))

need <- c("gj-exact-points-30seeds.csv", "gj-probing-30seeds.csv", "gj-feed-30seeds.csv",
          "gj-feedcost-30seeds.csv")
need_denv <- c("gj-denv-rateb-30seeds.csv", "gj-denv-probing-30seeds.csv",
               "gj-denv-points-30seeds.csv")
missing <- c(need[!file.exists(file.path(OUT, need))],
             need_denv[!file.exists(file.path(OUT_DENV, need_denv))])
if (length(missing)) stop("missing: ", paste(missing, collapse = ", "))

read_eip <- function(f, dir = OUT) read_complete(file.path(dir, f))

row <- function(panel, kind, series, x, pathogen, y, b, n = length(b))
  cbind(data.frame(panel = panel, kind = kind, series = series, x = x, pathogen = pathogen,
                   eip = EIP_DAYS, n_pairs = n, sum_manip = sum(y), sum_base = sum(b)),
        as.data.frame(t(pooled_fold(y, b))))

# sweep curves, each read from its own parameter set's folder
sweep_rows <- function(f, panel, series, pathogen, dir = OUT) {
  df <- read_eip(f, dir)
  base <- df[df$mult == min(df$mult), c("seed", "gj_vc")]; names(base)[2] <- "b"
  m <- merge(df, base, by = "seed"); m <- m[m$b > 0, ]
  do.call(rbind, lapply(sort(unique(m$mult)), function(mm) {
    z <- m[m$mult == mm, ]
    setNames(row(panel, "sweep", series, mm, pathogen, z$gj_vc, z$b),
             c("panel","kind","series","x","pathogen","eip",
               "n_pairs","sum_manip","sum_base","est","lo","hi"))
  }))
}

# exact scenario points, each from its own parameter set
points <- read_eip("gj-exact-points-30seeds.csv")
points_denv <- read_eip("gj-denv-points-30seeds.csv", OUT_DENV)
exact_rows <- function(scen, panel, pathogen, pts) {
  base_nm <- grep("^baseline", unique(pts$scenario), value = TRUE)[1]
  a <- pts[pts$scenario == scen, c("seed", "gj_vc")]
  b <- pts[pts$scenario == base_nm, c("seed", "gj_vc")]
  m <- merge(a, b, by = "seed"); m <- m[m$gj_vc.y > 0, ]
  setNames(row(panel, "exact", scen, NA_real_, pathogen, m$gj_vc.x, m$gj_vc.y),
           c("panel","kind","series","x","pathogen","eip","n_pairs","sum_manip","sum_base","est","lo","hi"))
}

# S6 arms
fc <- read_eip("gj-feedcost-30seeds.csv")
arm_rows <- function(arm, scen_name, anchor) {
  base <- anchor[, c("seed", "gj_vc")]; names(base)[2] <- "b"
  do.call(rbind, lapply(sort(unique(fc$feed_survival[fc$scenario == scen_name]), decreasing = TRUE),
    function(fs) {
      a <- fc[fc$scenario == scen_name & abs(fc$feed_survival - fs) < 1e-9, c("seed", "gj_vc")]
      m <- merge(a, base, by = "seed"); m <- m[m$b > 0, ]
      setNames(row("d", "arms", paste0("engorgement x1.92 + per-meal cost (", arm, ")"), 1 - fs,
                   "malaria", m$gj_vc, m$b),
               c("panel","kind","series","x","pathogen","eip","n_pairs","sum_manip","sum_base","est","lo","hi"))
    }))
}

res <- list(
  sweep_rows("gj-denv-rateb-30seeds.csv", "a", "landing rate (rate_b @ probeH 0.4)", "dengue", OUT_DENV),
  sweep_rows("gj-probing-30seeds.csv", "b", "probability of probing", "malaria"),
  sweep_rows("gj-denv-probing-30seeds.csv", "b", "probability of probing", "dengue", OUT_DENV),
  sweep_rows("gj-feed-30seeds.csv", "c", "probability of engorgement", "malaria"),
  exact_rows("S2_rate", "a", "dengue", points_denv),
  exact_rows("S3_probe", "b", "malaria", points),
  exact_rows("S4_probe", "b", "dengue", points_denv),
  exact_rows("S5_feed", "c", "malaria", points),
  arm_rows("boost", "feedcost", fc[fc$scenario == "baseline", ]),
  arm_rows("no boost", "feedcost_noboost",
           fc[fc$scenario == "feedcost_noboost" & abs(fc$feed_survival - 1) < 1e-9, ])
)
tab <- do.call(rbind, res)
tab <- tab[order(tab$panel, tab$kind, tab$series, tab$x, tab$pathogen), ]
out <- file.path(OUT, "gj-values-30seeds.csv")
write.csv(tab, out, row.names = FALSE, quote = FALSE)

cat(sprintf("=== plotted values at EIP %d d (estimand: sum VC(manip) / sum VC(base)) ===\n", EIP_DAYS))
for (i in seq_len(nrow(tab))) {
  xl <- if (is.na(tab$x[i])) "exact" else sprintf("%.2f", tab$x[i])
  cat(sprintf("  %-1s %-6s %-52s x = %-6s %6.3f [%6.3f, %6.3f]  n = %2d%s\n",
              tab$panel[i], tab$kind[i],
              substr(paste0(tab$series[i], " [", tab$pathogen[i], "]"), 1, 52), xl,
              tab$est[i], tab$lo[i], tab$hi[i], tab$n_pairs[i],
              if (tab$lo[i] > 1 || tab$hi[i] < 1) "  (CI excludes 1)" else ""))
}

cat("\n=== checks ===\n")
chk <- function(ok, msg) cat(if (isTRUE(ok)) "  PASS  " else "  FAIL  ", msg, "\n", sep = "")
sb <- tab[tab$kind == "sweep" & !is.na(tab$x) & tab$x == 1, ]
chk(nrow(sb) == 4 && all(abs(sb$est - 1) < 1e-12),
    "each sweep (both parameter sets) at x = 1.00 gives fold = 1.000 (by construction)")
nb1 <- tab[tab$kind == "arms" & grepl("no boost", tab$series) & abs(tab$x) < 1e-9, ]
chk(nrow(nb1) == 1 && abs(nb1$est - 1) < 1e-12, "no-boost arm at fs = 1.00 is exactly 1.000")
eng <- tab[tab$kind == "exact" & tab$series == "S5_feed", ]
bo <- tab[tab$kind == "arms" & grepl("\\(boost\\)$", tab$series) & abs(tab$x) < 1e-9, ]
chk(nrow(eng) == 1 && nrow(bo) == 1 && abs(eng$est - bo$est) < 1e-9,
    sprintf("engorgement is the same quantity in panels c and d (fold = %.3f in both)", eng$est[1]))
ms <- tab[tab$kind == "exact", ]
chk(nrow(ms) == 4 && all(ms$n_pairs == SEEDS_EXPECTED),
    sprintf("all four exact points pair %d seeds", SEEDS_EXPECTED))
cat("\nwrote", out, sprintf("(%d rows)\n", nrow(tab)))
