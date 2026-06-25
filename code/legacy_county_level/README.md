## Legacy County-Level Scripts

These scripts are retained for the earlier exploratory county-level pollution microbiome workflow. They are not part of the current SHRF/ZCTA analysis.

Do not run this folder for the current collaborator pipeline unless you specifically intend to reproduce the older county-level aggregate analyses.

Scripts:

1. `01_microbiome_cohort_export.R`
2. `02_pollution_microbe_correlation.R`
3. `03_plot_pollution_microbe_correlations.R`
4. `04_pollution_microbe_risk_models.R`
5. `05_hierarchical_microbe_phenotype_models.R`
6. `06_plot_hierarchical_findings.R`
7. `07_summarize_hierarchical_sensitivities.R`

Legacy requirements:

1. County-level exposure files configured with optional `exposome_path`.
2. County linkage through `hospitalization.county_code`.
3. Older culture windows and phenotype definitions documented in the scripts themselves.
