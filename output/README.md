## Output directory

Use this directory for generated local outputs. These files are produced by the scripts in `code/` and are intentionally separated from source code and documentation.

## Subdirectories

1. `final/`: aggregate CSV outputs from cohort export and exploratory correlation scripts.
2. `figures/`: generated PNG/PDF figures, including the pollution-microbe correlation tile plot.

## Expected Files

Current scripts may create:

1. `microbe_site_county_year_<site>_<stamp>.csv`
2. `microbe_organism_group_<site>_<stamp>.csv`
3. `microbe_organism_category_<site>_<stamp>.csv`
4. `microbe_qc_summary_<site>_<stamp>.csv`
5. `pollution_microbe_correlations_<stamp>.csv`
6. `pollution_microbe_correlation_tile_<stamp>.png`
7. `pollution_microbe_risk_models_<stamp>.csv`
8. `hierarchical_pollution_microbe_phenotype_models_<site>_<stamp>.csv`
9. `hierarchical_pollution_microbe_forest_<stamp>.png`
10. `hierarchical_pollution_microbe_heatmap_<stamp>.png`
11. `hierarchical_pollution_microbe_interactions_<stamp>.png`
12. `hierarchical_sensitivity_summary_<stamp>.csv`

## Governance

Generated site outputs should be treated as local analysis products. Do not commit patient-level data, CLIF source tables, site-specific configuration files, or unsuppressed small-cell aggregate outputs without explicit approval. For multi-site sharing, apply the release rules specified by the project team.

The repository can keep this README and placeholder folders under version control while leaving generated outputs local.
