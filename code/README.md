## Code directory

This directory contains the executable R workflow for the current CLIF pollution-microbiome project. The active pipeline is the severe hypoxemic respiratory failure (SHRF) ZCTA workflow. Scripts read site-specific paths from `config/config.json` through `utils/config.R` and `utils/clif_io.R`.

Older county-level exploratory scripts live in `code/legacy_county_level/`. They are retained for reference only and are not part of the current pipeline.

## Scripts

1. `08_count_severe_hypoxemic_rf.R`
   Counts adult ED-to-ICU hospitalizations with invasive mechanical ventilation and PaO2/FiO2 `<300` in the first 24 hours after ICU admission.

   Main output:
   - `severe_hypoxemic_rf_count_<site>_<stamp>.csv`

2. `09_count_pulmonary_cultures_in_shrf.R`
   Counts positive pulmonary cultures among SHRF hospitalizations. The current primary window is ICU admission through 48 hours after ICU admission.

   Main outputs:
   - `shrf_positive_pulmonary_culture_count_<site>_<stamp>.csv`
   - `shrf_positive_pulmonary_culture_top_organisms_<site>_<stamp>.csv`

3. `10_shrf_zcta_pollution_pulmonary_culture_models.R`
    Fits first-48h any-positive-pulmonary-culture mixed-effects logistic models against annualized PM2.5, annualized ozone, and annual NO2.

    Main outputs:
    - `shrf_zcta_pollution_any_pulmonary_culture_models_<site>_<stamp>.csv`
    - `shrf_zcta_pollution_any_pulmonary_culture_coverage_<site>_<stamp>.csv`

4. `11_shrf_zcta_pollution_organism_models.R`
    Fits first-48h organism-specific mixed-effects logistic models against annualized PM2.5, annualized ozone, and annual NO2.

    Main outputs:
    - `shrf_zcta_pollution_organism_models_<site>_<stamp>.csv`
    - `shrf_zcta_pollution_organism_model_counts_<site>_<stamp>.csv`

    Optional environment variables:
    - `MIN_ORGANISM_DETECTIONS`: minimum any-hospitalization pulmonary organism detections required for screening; default `25`.
    - `CULTURE_WINDOWS`: culture windows to model; default `first_48h`.

5. `12_plot_shrf_organism_forest.R`
    Creates a forest plot from the latest SHRF organism model output, with separate panels for PM2.5, ozone, and NO2.

    Optional environment variable:
    - `SHRF_ORGANISM_MODEL_PATH`: explicit organism-model CSV to plot.

## Run Order

From the repository root:

```bash
Rscript code/08_count_severe_hypoxemic_rf.R
Rscript code/09_count_pulmonary_cultures_in_shrf.R
Rscript code/10_shrf_zcta_pollution_pulmonary_culture_models.R
Rscript code/11_shrf_zcta_pollution_organism_models.R
Rscript code/12_plot_shrf_organism_forest.R
```

Sensitivity examples:

```bash
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

1. SHRF cohort uses adult ED-to-ICU hospitalizations with IMV and P/F `<300` in the first ICU 24 hours.
2. ED intubations are eligible.
3. SHRF/ZCTA models use cultures collected from ICU admission through 48 hours after ICU admission.
4. Primary respiratory specimens are `respiratory_tract` and `respiratory_tract_lower`.
5. Organism-specific models use `organism_category`.

## Development Notes

The current code intentionally keeps the exploratory outputs unsuppressed for local signal-finding. Before sharing outside an approved environment, add minimum-cell suppression or review aggregate release rules. The next major code step is to replace the respiratory support proxy with the fuller physiologic ARF phenotype using SpO2/FiO2, PaO2/FiO2, and PaCO2/pH logic.
