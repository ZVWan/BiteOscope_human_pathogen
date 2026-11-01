# Per-infected-female transmission readout

Counts, for every simulated mosquito, the infectious probes it delivers and the
human infections it causes, and turns them into the readout reported in the
manuscript: infectious bites per infected female (B, Cator et al. 2014).

The model itself is not modified. `metric-counter.R` rebinds two
`Mosquito_Female` methods after `MBITES/code/_common.R` (or
`MBITES/code/denv-pathogen.R`) has been sourced:

- `feedHost` records the time at which a mosquito first acquires the pathogen;
- `probeHost` counts probes by infectious mosquitoes, and the probes that
  actually infected the human.

Neither method draws random numbers, so a run with the counters returns the same
simulation as a run without them for the same seed and settings.

## Files

    metric-counter.R         the two rebound methods and metric_summary()
    plot-helpers.R           panel style and PDF device
    run_metric_points.R      exact scenarios (baseline, S2, S3, S4, S5), 30 seeds
    run_metric_all.R         sweeps 1-2.5x and the feed-cost analysis, 30 seeds
    analyse.R                fold change in B at day 60 and day 70
    plot_metric_panels.R     Fig. 4 panels
    results/                 stored runs and metric-values.csv

## Usage

    Rscript run_metric_points.R [n_workers]   # default 12
    Rscript run_metric_all.R    [n_workers]   # default 12
    Rscript analyse.R
    Rscript plot_metric_panels.R

The runners append every finished run to `results/metric-sweeps-30seeds.csv`
(or `metric-points-30seeds.csv`) and skip runs that are already stored, so an
interrupted job can be restarted. `run_metric_all.R` also stops when a file named
`STOP` is created in this folder.

`plot_metric_panels.R` writes `fig4{a,b,c,d}-metric.pdf` and `fig4-metric.pdf` to
`results/` (or to `METRIC_OUT`), together with `metric-values.csv`.

## What is computed

For each run, `metric_summary()` returns, at the 60- and 70-day cut-offs:

- `ninf_cut`: mosquitoes infected by that day;
- `probes_cut`: all probes by those mosquitoes that were infectious at the time,
  counted over the whole run;
- `tx_cut`: human infections those probes caused.

The readout is the pooled ratio over seeds, `sum(probes_cut) / sum(ninf_cut)`,
expressed as the fold change against the matched baseline on the same seeds
(ratio of ratios). Confidence intervals come from a paired percentile bootstrap
over the 30 seeds (B = 2000, seed 20261030 in the plotting script).

## Why the cut-off

Runs last 100 days and the incubation period is 11 days, so a female infected in
the final 11 days cannot become infectious before the run ends. Counting her as
having zero infectious bites would be administrative censoring, not biology.
Including only females infected by day 60 (or 70) gives every included female at
least 40 (or 30) days of follow-up, far beyond the incubation period; the
estimate is therefore conditional on infection by that day, and the two cut-offs
give the same values to two decimals.

## Results columns

- `probes_cut / ninf_cut`: infectious bites per infected female (B).
- `tx_cut`: human infections those probes caused, for reference.

## License

GNU General Public License v3.0, see `../LICENSE`.
