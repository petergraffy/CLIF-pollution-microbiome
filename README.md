# CLIF Pollution-Microbiome

## Acute MWAS development

**Second-site pilot:** [installation, synthetic test, local runner and aggregate return](docs/buddy_testing.md). Start with `Rscript code/35_buddy_smoke_test.R`, then configure the six core CLIF tables and run `Rscript code/34_run_buddy_site.R`.

The newest workflow is an admission-based time-stratified case-crossover screen
of acute daily PM2.5 and ozone exposure. It includes hospitalizations entering
the ICU within 24 hours, with respiratory cultures in the first 24/48/72 hours
of ICU entry. See [the acute MWAS run instructions](code/README.md#acute-pm25--ozone-mwas-current-development-workflow).
The [federated extension](code/README.md#federated-mwas-with-acs-ses-and-admission-diagnoses)
adds ACS ZCTA SES interactions, diagnosis-defined secondary outcomes, aggregate
site exports, and central meta-analysis. It now includes age, sex, recorded-race and clinical-severity interactions, a primary-cohort Table 1, and annual site summaries. This is the current development path;
its all-organism screen replaces the initial 50/100-count filters.
The older prior-year positive-culture screen below is retained for provenance.

## Overview

This project screens associations between acute residential ZCTA PM2.5/ozone exposure and respiratory culture-detected organisms in hospital admissions entering the ICU within 24 hours. The current workflow is a time-stratified case-crossover analysis with aggregate outputs for federated meta-analysis. The older prior-year pollution workflow is retained for provenance.

CLIF contains clinical microbiology culture and susceptibility data rather than sequencing-based microbiome assays. For that reason, this repository uses the phrase **respiratory microbial ecology** or **culture-detected organisms** rather than claiming to measure the full lung microbiome.

## CLIF Version

This project targets CLIF 2.1.

## Scientific Aims

1. Screen named respiratory organisms against daily PM2.5 and ozone, using seven-day exposure as primary and 3/14/28-day sensitivities.
2. Assess neighborhood SES modification, diagnosis-defined outcomes, culture selection and specimen/antibiotic sensitivities.
3. Validate inference and site portability with conditional count models, simulation calibration and federated aggregate pooling.

## Required CLIF tables and fields

Please refer to the [CLIF data dictionary](https://clif-icu.com/data-dictionary), [CLIF Tools](https://clif-icu.com/tools), [ETL Guide](https://clif-icu.com/etl-guide), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for more information on constructing the required tables and fields. 

The acute MWAS buddy pilot requires **patient, hospitalization, adt, microbiology_culture, medication_admin_intermittent and hospital_diagnosis**. See [the buddy guide](docs/buddy_testing.md) and run `Rscript code/33_site_preflight.R` for required fields. SOFA derivation is optional in the first pilot.

The following tables are used by the older prior-year PheWAS:
1. **patient**: `patient_id`, `sex_category`, `race_category`, `ethnicity_category`
2. **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `zipcode_five_digit`
3. **microbiology_culture**: `patient_id`, `hospitalization_id`, `organism_id`, `order_dttm`, `collect_dttm`, `result_dttm`, `fluid_name`, `fluid_category`, `method_name`, `method_category`, `organism_name`, `organism_category`, `organism_group`

Additional tables such as `adt`, `respiratory_support`, `labs`, `vitals`, `hospital_diagnosis`, and medication tables are used only by the legacy restricted-cohort workflows.

## Legacy exposure data

The current PheWAS workflow uses ZCTA-level exposure parquet files linked by `hospitalization.zipcode_five_digit` and the calendar year before admission. Set `zcta_exposure_dir` in `config/config.json` to a directory containing:

1. `air_pollution_zcta_pm25_monthly_2005_2023.parquet`
2. `air_pollution_zcta_no2_annual_2005_2025.parquet`

Monthly PM2.5 is annualized by ZIP/year in the landscape script. NO2 is already annual. Current models scale PM2.5 per 5 ug/m3 and NO2 per 10 ppb.

The older aggregate county-level workflow can still use `exposome_path` with county-year PM2.5/NO2 files, but that field is optional and not required for the current PheWAS.

## Current Cohort

The current PheWAS includes every hospitalization with at least one positive pulmonary culture:

1. Pulmonary specimens use `fluid_category` values `respiratory_tract` and `respiratory_tract_lower`.
2. Positive culture means `method_category == "culture"` and `organism_group` or `organism_category` is not `no_growth`.
3. The denominator for each organism model is hospitalizations with any positive pulmonary culture.
4. The outcome is whether a specific `organism_category` was present in that hospitalization.
5. Exposures are prior-year ZCTA PM2.5 and NO2 linked by five-digit ZIP code.

## Repository Layout

1. `code/`: active PheWAS scripts plus legacy workflow subfolders.
2. `config/`: site-specific runtime configuration template. Real `config.json` files are ignored.
3. `docs/`: project rationale, working definitions, and CLIF primer material.
4. `output/`: local generated aggregate outputs and figures. Site-derived output should not be committed unless explicitly approved.
5. `utils/`: shared config-loading utilities.

## Current PheWAS Workflow

Use this workflow for the current project. It does not require county-level exposure files.

1. Run `code/15_positive_lung_cultures_prior_year_pollution_landscape.R` to export all positive pulmonary cultures with prior-year ZCTA PM2.5 and NO2.
2. Run `code/17_plot_positive_lung_organism_phewas.R` to fit organism-specific models and create the main PheWAS plots.
3. Optionally run `code/16_plot_positive_lung_group_phewas.R` for organism-group-level companion plots.

The first-pass outputs are exploratory. They are intended to help assess signal and feasibility before adding site pooling, sensitivity cohorts, and richer clinical phenotypes.

## Legacy County-Level Workflow

Scripts `01`-`07` are retained in `code/legacy_county_level/` for earlier exploratory county-level analyses. The former SHRF/ZCTA and ED-to-ICU early-culture scripts are retained in `code/legacy_shrf_zcta/`. Do not run these folders for the current PheWAS unless you specifically intend to reproduce older workflows.

## Outputs

The current PheWAS scripts save these files in [`output/final`](output/README.md) and [`output/figures`](output/README.md):

1. `positive_lung_cultures_prior_year_pollution_<site>_<stamp>.csv`
2. `positive_lung_cultures_prior_year_pollution_coverage_<site>_<stamp>.csv`
3. `positive_lung_cultures_prior_year_pollution_organism_summary_<site>_<stamp>.csv`
4. `positive_lung_cultures_prior_year_pollution_year_summary_<site>_<stamp>.csv`
5. `positive_lung_culture_organism_prior_year_pollution_models_<stamp>.csv`
6. `positive_lung_culture_organism_phewas_pm25_prior_year_<stamp>.png`
7. `positive_lung_culture_organism_phewas_no2_prior_year_<stamp>.png`

The legacy county-level scripts produce:

1. `microbe_site_county_year_<site>_<stamp>.csv`
2. `microbe_organism_group_<site>_<stamp>.csv`
3. `microbe_organism_category_<site>_<stamp>.csv`
4. `pollution_microbe_correlations_<stamp>.csv`
5. `pollution_microbe_risk_models_<stamp>.csv`
6. `hierarchical_pollution_microbe_phenotype_models_<site>_<stamp>.csv`

See [`docs/project_spec.md`](docs/project_spec.md) for the full working analysis plan.

## Legacy PheWAS run instructions

### 1. Update `config/config.json`

Copy [`config/config_template.json`](config/config_template.json) to `config/config.json` and update the site name, CLIF table path, file type, repository path, and ZCTA exposure directory. Follow the notes in [config/README.md](config/README.md).

### 2. Set up the project environment

This project uses R and `renv`. From the repository root:

```r
renv::restore()
```

If you already have the required packages installed, the scripts can be run directly with `Rscript`.

### 3. Run code

```bash
Rscript code/15_positive_lung_cultures_prior_year_pollution_landscape.R
Rscript code/17_plot_positive_lung_organism_phewas.R
```

Sensitivity examples:

```bash
MIN_ORGANISM_DETECTIONS=10 Rscript code/17_plot_positive_lung_organism_phewas.R
MIN_GROUP_DETECTIONS=10 Rscript code/16_plot_positive_lung_group_phewas.R
```

Detailed workflow instructions are provided in the [code directory](code/README.md).

## Data Governance

Do not commit patient-level CLIF tables, site configs, or unsuppressed site-derived outputs. The local exploratory output files can remain in `output/`, but they should be reviewed for sharing rules before being pushed or distributed.

## Next steps

Run the synthetic smoke test, then execute the frozen buddy protocol at another CLIF site. Review source mapping, coverage, sparse fits and simulation error rates locally before releasing institutionally approved aggregate outputs. Development calibration has identified excess null rejection in some scenarios; inferential validity remains an explicit research task.
