# CLIF Pollution-Microbiome

## Overview

This project studies whether ambient air pollution is associated with geographic variation in respiratory microbial ecology among adult ICU patients in CLIF. The core idea is that exposures such as PM2.5 and NO2 may shape which respiratory organisms are recovered from clinical cultures, and those culture-detected microbial profiles may relate to acute respiratory failure (ARF), pneumonia, sepsis, and respiratory support severity.

CLIF contains clinical microbiology culture and susceptibility data rather than sequencing-based microbiome assays. For that reason, this repository uses the phrase **respiratory microbial ecology** or **culture-detected organisms** rather than claiming to measure the full lung microbiome.

## CLIF Version

This project targets CLIF 2.1.

## Scientific Aims

1. Estimate whether county-year PM2.5 and NO2 exposures are associated with respiratory culture organism composition among adult ICU hospitalizations.
2. Test whether pollution-associated organism profiles are associated with ARF vulnerability, ARF severity, and respiratory support escalation.
3. Explore pneumonia- and sepsis-enriched subcohorts for distinct microbial-respiratory failure phenotypes.

## Required CLIF tables and fields

Please refer to the [CLIF data dictionary](https://clif-icu.com/data-dictionary), [CLIF Tools](https://clif-icu.com/tools), [ETL Guide](https://clif-icu.com/etl-guide), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for more information on constructing the required tables and fields. 

The following tables are required:
1. **patient**: `patient_id`, `sex_category`, `race_category`, `ethnicity_category`, `death_dttm`
2. **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `discharge_category`, `county_code`
3. **adt**: `hospitalization_id`, `hospital_id`, `in_dttm`, `out_dttm`, `location_category`
4. **microbiology_culture**: `patient_id`, `hospitalization_id`, `organism_id`, `collect_dttm`, `fluid_category`, `method_category`, `organism_category`, `organism_group`
5. **microbiology_susceptibility**: `organism_id`, `antimicrobial_category`, `susceptibility_category`; optional for exploratory culture-only analyses
6. **respiratory_support**: `hospitalization_id`, `recorded_dttm`, `device_category`

Recommended tables for the full analysis:
1. **vitals**: SpO2 for physiologic ARF.
2. **labs**: PaO2, PaCO2, and pH for hypoxemic/hypercapnic ARF.
3. **hospital_diagnosis**: pneumonia and sepsis subcohorts.
4. **medication_admin_continuous**: vasopressors and sepsis severity covariates.

## Exposure Data

The current workflow expects county-year exposure files with county FIPS and year columns. It can read these filenames when present in `config$exposome_path`:

1. `conus_county_pm25_2005_2024.csv` or `pm25_county_year.csv`
2. `conus_county_no2_2005_2024.csv` or `no2_county_year.csv`

The default linkage is `hospitalization.county_code` by admission year. Future analyses can add monthly exposure windows, lagged annual exposures, weather, SVI, and other county-level covariates.

## Cohort identification

Primary cohort:
1. Adult hospitalizations, age >= 18.
2. ICU admission identified through `adt.location_category == "icu"`.
3. ICU length of stay >= 24 hours.
4. Valid CONUS county FIPS in `hospitalization.county_code`.
5. Admission/ICU year covered by county-year PM2.5 and NO2 files.

Primary respiratory culture window:
1. 48 hours before through 72 hours after first ICU admission.
2. Respiratory specimens: `respiratory_tract` and `respiratory_tract_lower`.
3. Sensitivity specimens: add upper airway/oropharynx and pleural fluid.

## Repository Layout

1. `code/`: R scripts for cohort export, pollution-microbe correlations, and exploratory plots.
2. `config/`: site-specific runtime configuration template. Real `config.json` files are ignored.
3. `docs/`: project rationale, working definitions, and CLIF primer material.
4. `output/`: local generated aggregate outputs and figures. Site-derived output should not be committed unless explicitly approved.
5. `utils/`: shared config-loading utilities.

## Current Workflow

1. Run `code/01_microbiome_cohort_export.R` to build adult ICU respiratory culture aggregates and link county-year pollution.
2. Run `code/02_pollution_microbe_correlation.R` to compute crude county-year correlations between PM2.5/NO2 and organism-category prevalence.
3. Run `code/03_plot_pollution_microbe_correlations.R` to create a correlation tile plot.
4. Run `code/04_pollution_microbe_risk_models.R` to fit aggregate binomial models estimating organism detection odds per IQR increase in PM2.5 or NO2.
5. Run `code/05_hierarchical_microbe_phenotype_models.R` to fit patient-level mixed-effects logistic models with county random intercepts and pneumonia, sepsis, or severe respiratory support phenotype interactions.
6. Run `code/06_plot_hierarchical_findings.R` and `code/07_summarize_hierarchical_sensitivities.R` to create figures and compare sensitivity runs.

The first-pass outputs are exploratory. They are intended to help assess signal and feasibility before adding adjusted models, ARF physiology, pneumonia/sepsis definitions, and multi-site pooling.

## Outputs

The export script saves aggregate files in [`output/final`](output/README.md):

1. `microbe_site_county_year_<site>_<stamp>.csv`
2. `microbe_organism_group_<site>_<stamp>.csv`
3. `microbe_organism_category_<site>_<stamp>.csv`
4. `microbe_qc_summary_<site>_<stamp>.csv`

The correlation and plotting scripts produce:

1. `pollution_microbe_correlations_<stamp>.csv`
2. `pollution_microbe_correlation_tile_<stamp>.png`
3. `pollution_microbe_risk_models_<stamp>.csv`
4. `hierarchical_pollution_microbe_phenotype_models_<site>_<stamp>.csv`
5. `hierarchical_pollution_microbe_forest_<stamp>.png`
6. `hierarchical_pollution_microbe_heatmap_<stamp>.png`
7. `hierarchical_sensitivity_summary_<stamp>.csv`

See [`docs/project_spec.md`](docs/project_spec.md) for the full working analysis plan.

## Detailed Instructions for running the project

### 1. Update `config/config.json`

Copy [`config/config_template.json`](config/config_template.json) to `config/config.json` and update the site name, CLIF table path, file type, repository path, and exposome path. Follow the notes in [config/README.md](config/README.md).

### 2. Set up the project environment

This project uses R and `renv`. From the repository root:

```r
renv::restore()
```

If you already have the required packages installed, the scripts can be run directly with `Rscript`.

### 3. Run code

```bash
Rscript code/01_microbiome_cohort_export.R
Rscript code/02_pollution_microbe_correlation.R
Rscript code/03_plot_pollution_microbe_correlations.R
Rscript code/04_pollution_microbe_risk_models.R
Rscript code/05_hierarchical_microbe_phenotype_models.R
Rscript code/06_plot_hierarchical_findings.R
Rscript code/07_summarize_hierarchical_sensitivities.R
```

Sensitivity examples:

```bash
ANALYSIS_LABEL=no_cook_county EXCLUDE_COUNTY_FIPS=17031 Rscript code/05_hierarchical_microbe_phenotype_models.R
ANALYSIS_LABEL=lower_resp_only RESP_FLUID_MODE=lower_only MIN_ORGANISM_DETECTIONS=25 Rscript code/05_hierarchical_microbe_phenotype_models.R
```

Detailed workflow instructions are provided in the [code directory](code/README.md).

## Data Governance

Do not commit patient-level CLIF tables, site configs, or unsuppressed site-derived outputs. The local exploratory output files can remain in `output/`, but they should be reviewed for sharing rules before being pushed or distributed.

## Next Steps

1. Port the physiologic ARF phenotype from prior CLIF ARF pollution work.
2. Add pneumonia/sepsis subcohort definitions using `hospital_diagnosis` and optional medication/lab criteria.
3. Add sensitivity analyses for lower respiratory specimens only, positive cultures only, and Cook County versus non-Cook County catchments.
4. Extend the hierarchical models to multi-site pooled analysis with site random effects.
