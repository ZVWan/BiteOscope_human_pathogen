# Fold change in the per-infected-female readout against the baseline on the same
# seeds: infectious probes per infected female (females infected by day `cut`).
# Ratios are pooled over seeds, CI = paired percentile bootstrap over seeds (B = 2000).

# Paths are resolved from this script, so the package runs from any checkout.
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
d <- read.csv(file.path(HERE, "results", "metric-points-30seeds.csv"))

fold <- function(m, b, num, den, B = 2000) {
  m <- m[match(b$seed, m$seed), ]; ok <- !is.na(m$seed); m <- m[ok, ]; b <- b[ok, ]
  r <- function(i) (sum(m[i, num]) / sum(m[i, den])) / (sum(b[i, num]) / sum(b[i, den]))
  set.seed(20261030)
  bs <- replicate(B, r(sample(nrow(m), replace = TRUE)))
  sprintf("%.2f [%.2f-%.2f]", r(seq_len(nrow(m))), quantile(bs, .025), quantile(bs, .975))
}

for (cut in c(70, 60)) {
  cat(sprintf("\n== females infected by day %d ==\n", cut))
  res <- NULL
  for (set in unique(d$set)) {
    z <- d[d$set == set, ]; bn <- grep("^baseline", unique(z$scenario), value = TRUE)
    b <- z[z$scenario == bn, ]
    for (nm in setdiff(unique(z$scenario), bn)) {
      m <- z[z$scenario == nm, ]
      res <- rbind(res, data.frame(settings = set, scenario = nm, n = nrow(m),
        B = fold(m, b, paste0("probes_", cut), paste0("ninf_", cut))))
    }
  }
  print(res, row.names = FALSE)
}
