# MBITES model code and transmission readouts

Code for the transmission-model part of the manuscript on pathogen-induced changes
in mosquito behaviour. The model and the analysis of its transmission output live
in one tree, so a single checkout is enough to reproduce the reported numbers.

## Layout

    MBITES/
      code/
        _common.R                    settings, run_sim(), EIP pinning
        denv-pathogen.R              dengue pathogen settings and run_denv()
        manuscript/
          _senescence-common.R       model configuration and landscape
          manipulation-enhanced-probing.R
                                     infection-dependent probing, engorgement and host-seeking

    metric/
      metric-counter.R               counts infectious probes and human infections per mosquito
      plot-helpers.R                 plotting helpers
      run_metric_points.R            exact scenarios, 30 seeds
      run_metric_all.R               full sweeps and the feed-cost analysis
      analyse.R                      fold change in the readout, day-60 and day-70
      plot_metric_panels.R           Figures used in the manuscript
      results/                       stored runs and the plotted values

## Requirements

R 4.x, the MBITES R package (not on CRAN) and three CRAN packages:

    install.packages(c("futile.logger", "jsonlite", "ggplot2"))
    remotes::install_github("dd-harp/MBITES", ref = "v1.0", subdir = "MBITES")

The runners default to 12 workers; pass another number as the only argument,
e.g. `Rscript run_metric_all.R 20`.

## Running

    Rscript metric/run_metric_points.R
    Rscript metric/run_metric_all.R
    Rscript metric/analyse.R
    Rscript metric/plot_metric_panels.R

Every script locates itself, so it can be called from any directory. The metric
scripts use the model tree in `MBITES/`; set `MBITES_MODEL_DIR` to point at a
different copy.

`metric/results/` holds the per-run counters and `metric-values.csv`, the only
fold-change table.

## Readouts

**Infectious bites per infected female (B).** The readout used in the manuscript,
following Cator et al. (2014): the average number of infectious bites an infected
female delivers in her lifetime, the product of the probability of surviving to
the transmissible stage of infection and the average number of infectious bites
per mosquito reaching that stage. It is estimated from the simulated biting
histories as the total number of infectious probes divided by the number of
infected females, so females that die before becoming infectious contribute zero.
A probe by an infectious mosquito counts as an infectious bite because
transmission in MBITES occurs when the mosquito probes the host and does not
require a completed blood meal.

## Model settings

- 1000 female mosquitoes, 200 humans, 10% initial human prevalence, 100 days,
  incubation period 11 days, senescence on (sns_a = 0.25, sns_b = 60, onset
  about 16 days).
- Seeds 101, 202, ..., 3030 (30 replicates).
- Dengue runs: b 0.55, c 0.50, human incubation 6 days, human infectious period
  5 days; the incubation period stays at 11 days. Values and sources are in
  `MBITES/code/denv-pathogen.R`.
- `feed_survival` is added in `manipulation-enhanced-probing.R` as the per-meal
  mortality of infectious mosquitoes and is used only for the feed-cost analysis.

## License

GNU General Public License v3.0, see `LICENSE`.

The model is built on the MBITES R package (Wu et al., PLoS Comput Biol 2020,
doi:10.1371/journal.pcbi.1007446). The snapshot in `MBITES/code/manuscript/` is
derived from that package and is included for reproducibility with attribution;
copyright in that code remains with its authors. The GPL-3.0 grant covers the
code written for this study.
