## Legacy SHRF and Early-Culture ZCTA Scripts

These scripts are retained for the earlier severe hypoxemic respiratory failure
and ED-to-ICU early respiratory culture analyses. They are not part of the
current positive-lung-culture PheWAS workflow.

Do not run this folder for buddy testing the current PheWAS unless you
specifically intend to reproduce the older cohort-restricted analyses.

Scripts:

1. `08_count_severe_hypoxemic_rf.R`
2. `09_count_pulmonary_cultures_in_shrf.R`
3. `10_shrf_zcta_pollution_pulmonary_culture_models.R`
4. `11_shrf_zcta_pollution_organism_models.R`
5. `12_plot_shrf_organism_forest.R`
6. `13_count_ed_icu_early_resp_culture_cohort.R`
7. `14_ed_icu_early_resp_culture_pollution_organism_models.R`

Legacy requirements:

1. CLIF 2.1 tables configured through `config/config.json` or environment
   overrides.
2. ZCTA exposure files configured with `zcta_exposure_dir` or
   `ZCTA_EXPOSURE_DIR`.
3. Additional CLIF tables for respiratory support, labs, ADT, and medication
   logic depending on the script.
