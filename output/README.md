## Output Directory

Use this directory for generated local outputs. These files are produced by the
scripts in `code/` and are intentionally separated from source code and
documentation.

Generated site outputs are ignored by default. Do not commit patient-level data,
site-specific CLIF extracts, or unsuppressed site-derived outputs unless the
project team explicitly approves release.

## Subdirectories

1. `final/`: CSV exports, summaries, and model results.
2. `figures/`: generated PNG/PDF figures.

## Current PheWAS Outputs

The active positive-lung-culture PheWAS workflow may create:

1. `positive_lung_cultures_prior_year_pollution_<site>_<stamp>.csv`
2. `positive_lung_cultures_prior_year_pollution_coverage_<site>_<stamp>.csv`
3. `positive_lung_cultures_prior_year_pollution_organism_summary_<site>_<stamp>.csv`
4. `positive_lung_cultures_prior_year_pollution_year_summary_<site>_<stamp>.csv`
5. `positive_lung_culture_organism_prior_year_pollution_models_<stamp>.csv`
6. `positive_lung_culture_organism_phewas_site_summary_<stamp>.csv`
7. `positive_lung_culture_organism_phewas_modeled_organisms_<stamp>.csv`
8. `positive_lung_culture_organism_phewas_pm25_prior_year_<stamp>.png`
9. `positive_lung_culture_organism_phewas_no2_prior_year_<stamp>.png`

Optional group-level companion outputs:

1. `positive_lung_culture_group_prior_year_pollution_models_<stamp>.csv`
2. `positive_lung_culture_group_phewas_site_summary_<stamp>.csv`
3. `positive_lung_culture_group_phewas_modeled_groups_<stamp>.csv`
4. `positive_lung_culture_group_phewas_pm25_prior_year_<stamp>.png`
5. `positive_lung_culture_group_phewas_no2_prior_year_<stamp>.png`

## Legacy Outputs

Older county-level and restricted-cohort workflows may create files beginning
with `microbe_`, `pollution_microbe_`, `hierarchical_`, `shrf_`, or
`ed_icu_early_resp_`. Those are retained locally for provenance but are not part
of the current buddy-test workflow.
