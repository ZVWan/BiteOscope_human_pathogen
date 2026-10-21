# MBITES model code: infection-dependent behavioural manipulation

Code for the transmission-model part of the manuscript on pathogen-induced
changes in mosquito behaviour.

## Layout

    code/                       
      _senescence-common.R         model settings, landscape and Garrett-Jones VC readout
      manipulation-enhanced-probing.R
                                   rebinds Mosquito_Female methods so that
                                   probing, engorgement and host-seeking rate
                                   depend on the mosquito's infection status

    30seeds_sweeps/              
      run_sweep30.R                three Garrett-Jones sweeps, 30 seeds per point
      run_exact_scenario.R         exact points S2/S3/S4/S5 plus shared baseline
      run_s6_feedcost.R            S6 trade-off: engorgement boost and survival cost
      run_s6_feedcost_noboost.R    same sweep without the boost
      plot_pooled_sweeps30.R       Fig. 4a-c panels
      plot_feedcost_arms30.R       Fig. 4d panel
      summary_pooled_values.R      plotted values in one table

## Running

Each script locates itself, so it can be called from anywhere. Run
the scripts in 30seeds_sweeps/ rather than sourcing files in code/ directly.

    Rscript 30seeds_sweeps/run_sweep30.R 20        
    Rscript 30seeds_sweeps/run_exact_scenario.R 20
    Rscript 30seeds_sweeps/run_s6_feedcost.R 20
    Rscript 30seeds_sweeps/run_s6_feedcost_noboost.R 20
    Rscript 30seeds_sweeps/plot_pooled_sweeps30.R
    Rscript 30seeds_sweeps/plot_feedcost_arms30.R
    Rscript 30seeds_sweeps/summary_pooled_values.R


## Settings

  - 1000 female mosquitoes, 200 humans, 10% initial human prevalence, 100 days,
    EIP 14 days, senescence on (sns_a = 0.25, sns_b = 60, onset about 16 days).
  - Seeds 101, 202, ..., 3030 (30 replicates).
  - Reported effect: fold change in Garrett-Jones vectorial capacity against a
    matched uninfected baseline, pooled over the paired seeds
    (sum(VC manipulated) / sum(VC baseline)), with a paired percentile bootstrap
    95% CI (B = 2000).
  - feed_survival is added in manipulation-enhanced-probing.R to represent the per-meal mortality of
    infectious mosquitoes and is only used for the Fig. 4d analysis.
