# MBITES model code: infection-conditional behavioral manipulation

Code and results for the transmission-model part of the manuscript on
pathogen-induced behavioral changes in mosquitoes.

## Layout

    code/                        sourced by every script below
      _common.R                    settings, the run function, the vectorial
                                   capacity readout and the CSV storage
      denv-pathogen.R              dengue pathogen settings for the denv_sweeps scripts
      manuscript/                  the model code, copied from the repository:
        _senescence-common.R         model configuration, landscape and the
                                     Garrett-Jones VC readout
        manipulation-enhanced-probing.R
                                     rebinds the Mosquito_Female methods so that
                                     probing, engorgement and the host-seeking
                                     rate depend on the infection status of the
                                     individual mosquito

    30seeds_sweeps/              scripts producing the Fig. 4 values
      run_sweeps.R                 three Garrett-Jones sweeps, 30 seeds per point
      run_exact_points.R           all exact points S2/S3/S4/S5 + shared baseline
      run_feedcost.R               Trade-off (S6): engorgement boost + survival cost
      plot_panels.R                Fig. 4a-c panels
      plot_feedcost.R              Fig. 4d panel
      summary_values.R             every plotted value in one table

    denv_sweeps/                 the same scenarios under the dengue pathogen settings
      run_denv_sweeps.R            probing and landing-rate sweeps
      run_denv_points.R            exact points S2 and S4 + shared baseline
      plot_denv_panels.R           the two dengue panels
      plot_probing_malaria_dengue.R  probing under both pathogen settings in one figure

## Running

Each script locates itself, so it can be called from anywhere. The sweep scripts pass the code folder to their workers through
MBITES_COMMON_DIR, so source the scripts rather than the files in code/.

    Rscript 30seeds_sweeps/run_sweeps.R 20   
    Rscript 30seeds_sweeps/run_exact_points.R 20
    Rscript 30seeds_sweeps/run_feedcost.R 20
    Rscript denv_sweeps/run_denv_sweeps.R 20
    Rscript denv_sweeps/run_denv_points.R 20
    Rscript 30seeds_sweeps/plot_panels.R
    Rscript 30seeds_sweeps/plot_feedcost.R
    Rscript denv_sweeps/plot_denv_panels.R
    Rscript 30seeds_sweeps/plot_probing_malaria_dengue.R
    Rscript 30seeds_sweeps/summary_values.R

## Settings

  - 1000 female mosquitoes, 200 humans, 10% initial human prevalence, 100 days,
    EIP 11 days, senescence on (sns_a = 0.25, sns_b = 60, onset about 16 days).
  - Seeds 101, 202, ..., 3030 (30 replicates).
  - Reported effect: fold change in Garrett-Jones vectorial capacity against a
    matched uninfected baseline, pooled over the paired seeds
    (sum(VC manipulated) / sum(VC baseline)), with a paired percentile bootstrap
    95% CI (B = 2000).
  - feed_survival is not a native MBITES parameter. It is added in
    manipulation-enhanced-probing.R to represent the per-meal mortality of
    infectious mosquitoes and is only used for the Fig. 4d analysis.
  - Dengue runs (denv_sweeps/): b 0.55, c 0.50, human incubation 6 days and a
    human infectious period of 5 days, so an infection clears on day 11. The
    values and their sources are in code/denv-pathogen.R.
