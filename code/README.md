## Code directory

This directory contains the executable R workflow for the CLIF pollution-microbiome project. The implementation includes the original aggregate county-year workflow and a severe hypoxemic respiratory failure (SHRF) ZCTA workflow. Scripts read site-specific paths from `config/config.json` through `utils/config.R` and `utils/clif_io.R`.

## Scripts

1. `01_microbiome_cohort_export.R`
   Builds the adult ICU cohort, identifies early respiratory cultures, summarizes organism groups and organism categories, links optional susceptibility data, links county-year exposures, and writes aggregate CSVs.

   Main outputs:
   - `microbe_site_county_year_<site>_<stamp>.csv`
   - `microbe_organism_group_<site>_<stamp>.csv`
   - `microbe_organism_category_<site>_<stamp>.csv`
   - `microbe_qc_summary_<site>_<stamp>.csv`

2. `02_pollution_microbe_correlation.R`
   Reads the latest export files, completes the organism-by-county-year grid so absent organisms are treated as zero counts, and computes Pearson and Spearman correlations for PM2.5 and NO2.

   Main output:
   - `pollution_microbe_correlations_<stamp>.csv`

3. `03_plot_pollution_microbe_correlations.R`
   Reads the latest correlation output and creates a tile plot of organism-category correlations with pollution.

   Main outputs:
   - `pollution_microbe_correlation_tile_<stamp>.png`
   - `pollution_microbe_correlation_tile_<stamp>.pdf`

4. `04_pollution_microbe_risk_models.R`
   Fits aggregate quasi-binomial models for each organism category. Models estimate the odds ratio for organism detection per IQR increase in PM2.5 or NO2, adjusting for admission year. The script evaluates both the respiratory-cultured denominator and the positive-respiratory-culture denominator.

   Main output:
   - `pollution_microbe_risk_models_<stamp>.csv`

5. `05_hierarchical_microbe_phenotype_models.R`
   Rebuilds the patient-level respiratory culture cohort locally and fits mixed-effects logistic models for top organisms. Each model includes a county random intercept, year, age band, sex, pollutant exposure, phenotype, and pollutant-by-phenotype interaction.

   Main output:
   - `hierarchical_pollution_microbe_phenotype_models_<site>_<stamp>.csv`

   Optional environment variables:
   - `ANALYSIS_LABEL`: suffix for the output filename.
   - `EXCLUDE_COUNTY_FIPS`: comma-separated county FIPS values to exclude.
   - `MIN_ORGANISM_DETECTIONS`: minimum organism detections required for modeling.
   - `RESP_FLUID_MODE=lower_only`: restricts to `respiratory_tract_lower` specimens.

6. `06_plot_hierarchical_findings.R`
   Creates forest, heatmap, and interaction figures from the main hierarchical model output.

7. `07_summarize_hierarchical_sensitivities.R`
   Summarizes all hierarchical model outputs in `output/final` and writes sensitivity comparison tables.

8. `08_count_severe_hypoxemic_rf.R`
   Counts adult ED-to-ICU hospitalizations with invasive mechanical ventilation and PaO2/FiO2 `<300` in the first 24 hours after ICU admission.

   Main output:
   - `severe_hypoxemic_rf_count_<site>_<stamp>.csv`

9. `09_count_pulmonary_cultures_in_shrf.R`
   Counts positive pulmonary cultures among SHRF hospitalizations. The current primary window is ICU admission through 48 hours after ICU admission.

   Main outputs:
   - `shrf_positive_pulmonary_culture_count_<site>_<stamp>.csv`
   - `shrf_positive_pulmonary_culture_top_organisms_<site>_<stamp>.csv`

10. `10_shrf_zcta_pollution_pulmonary_culture_models.R`
    Fits first-48h any-positive-pulmonary-culture mixed-effects logistic models against annualized PM2.5, annualized ozone, and annual NO2.

    Main outputs:
    - `shrf_zcta_pollution_any_pulmonary_culture_models_<site>_<stamp>.csv`
    - `shrf_zcta_pollution_any_pulmonary_culture_coverage_<site>_<stamp>.csv`

11. `11_shrf_zcta_pollution_organism_models.R`
    Fits first-48h organism-specific mixed-effects logistic models against annualized PM2.5, annualized ozone, and annual NO2.

    Main outputs:
    - `shrf_zcta_pollution_organism_models_<site>_<stamp>.csv`
    - `shrf_zcta_pollution_organism_model_counts_<site>_<stamp>.csv`

    Optional environment variables:
    - `MIN_ORGANISM_DETECTIONS`: minimum any-hospitalization pulmonary organism detections required for screening; default `25`.
    - `CULTURE_WINDOWS`: culture windows to model; default `first_48h`.

12. `12_plot_shrf_organism_forest.R`
    Creates a forest plot from the latest SHRF organism model output, with separate panels for PM2.5, ozone, and NO2.

    Optional environment variable:
    - `SHRF_ORGANISM_MODEL_PATH`: explicit organism-model CSV to plot.

## Run Order

From the repository root:

```bash
Rscript code/01_microbiome_cohort_export.R
Rscript code/02_pollution_microbe_correlation.R
Rscript code/03_plot_pollution_microbe_correlations.R
Rscript code/04_pollution_microbe_risk_models.R
Rscript code/05_hierarchical_microbe_phenotype_models.R
Rscript code/06_plot_hierarchical_findings.R
Rscript code/07_summarize_hierarchical_sensitivities.R
Rscript code/08_count_severe_hypoxemic_rf.R
Rscript code/09_count_pulmonary_cultures_in_shrf.R
Rscript code/10_shrf_zcta_pollution_pulmonary_culture_models.R
Rscript code/11_shrf_zcta_pollution_organism_models.R
Rscript code/12_plot_shrf_organism_forest.R
```

Sensitivity examples:

```bash
ANALYSIS_LABEL=no_cook_county EXCLUDE_COUNTY_FIPS=17031 Rscript code/05_hierarchical_microbe_phenotype_models.R
ANALYSIS_LABEL=no_cook_county_min25 EXCLUDE_COUNTY_FIPS=17031 MIN_ORGANISM_DETECTIONS=25 Rscript code/05_hierarchical_microbe_phenotype_models.R
ANALYSIS_LABEL=lower_resp_only RESP_FLUID_MODE=lower_only MIN_ORGANISM_DETECTIONS=25 Rscript code/05_hierarchical_microbe_phenotype_models.R
MIN_ORGANISM_DETECTIONS=25 Rscript code/11_shrf_zcta_pollution_organism_models.R
SHRF_ORGANISM_MODEL_PATH=output/final/shrf_zcta_pollution_organism_models_YOUR_SITE_YYYYMMDD_HHMMSS.csv Rscript code/12_plot_shrf_organism_forest.R
```

## Configuration

Copy `config/config_template.json` to `config/config.json` and set:

1. `site_name`: short label for output filenames.
2. `tables_path`: directory containing CLIF tables.
3. `file_type`: `parquet`, `csv`, `fst`, or `auto`.
4. `zcta_exposure_dir`: directory containing the ZCTA PM2.5, ozone, and NO2 parquet release files.

The SHRF scripts locate CLIF tables recursively under `tables_path` and accept filenames with or without the `clif_` prefix.

## Current Analytic Definitions

1. Base cohort: adult ICU hospitalizations with ICU length of stay at least 24 hours and a valid CONUS county FIPS.
2. Culture window: 48 hours before through 72 hours after first ICU admission.
3. Primary respiratory specimens: `respiratory_tract` and `respiratory_tract_lower`.
4. Organism group analysis uses `organism_group`.
5. Genus/species analysis uses `organism_category`.
6. Severity proxy uses early HFNC, NIV/NIPPV/CPAP, or IMV in `respiratory_support`.
7. SHRF cohort uses ED-to-ICU hospitalizations with IMV and P/F `<300` in the first ICU 24 hours.
8. SHRF/ZCTA models use cultures collected from ICU admission through 48 hours after ICU admission.

## Development Notes

The current code intentionally keeps the exploratory outputs unsuppressed for local signal-finding. Before sharing outside an approved environment, add minimum-cell suppression or review aggregate release rules. The next major code step is to replace the respiratory support proxy with the fuller physiologic ARF phenotype using SpO2/FiO2, PaO2/FiO2, and PaCO2/pH logic.
