# MBITES model code: infection-dependent behavioural manipulation

Model setup for the transmission part of the manuscript: configuration, landscape
and the two methods that make probing, engorgement and host-seeking depend on the
mosquito's infection status. The runs, the readout and the figures are in
`../metric/`.

## Requirements

    install.packages(c("futile.logger", "jsonlite", "ggplot2"))
    remotes::install_github("dd-harp/MBITES", ref = "v1.0", subdir = "MBITES")

## Layout

    code/
      _common.R                    settings, run_sim(), EIP pinning
      denv-pathogen.R              dengue pathogen settings and run_denv()
      manuscript/
        _senescence-common.R       model configuration and landscape
        manipulation-enhanced-probing.R
                                   infection-dependent probing, engorgement and host-seeking

## Settings

- 1000 female mosquitoes, 200 humans, 10% initial human prevalence, 100 days,
  EIP 11 days, senescence on (sns_a = 0.25, sns_b = 60, onset about 16 days).
- Seeds 101, 202, ..., 3030 (30 replicates).
- `feed_survival` is added in `manipulation-enhanced-probing.R` as the per-meal
  mortality of infectious mosquitoes and is used only for the feed-cost analysis.
- Dengue runs: b 0.55, c 0.50, human incubation 6 days, human infectious period
  5 days; the incubation period stays at 11 days. Values and sources are in
  `code/denv-pathogen.R`.

## Use

    SCRIPT_DIR <- "MBITES"          # the folder that holds code/
    source(file.path(SCRIPT_DIR, "code", "_common.R"))
    run_sim(101)

`_common.R` stops a run if the two files in `code/manuscript/` change: the MD5
values it checks are the revision the stored results came from.
