## Code Directory

This directory contains the executable R workflow for the current CLIF
pollution-microbiome project. The active workflow is now the positive lung
culture organism-wide pollution association screen, or PheWAS-style workflow.

Scripts read site-specific paths from `config/config.json` through
`utils/config.R` and `utils/clif_io.R`. They can also be buddy tested with
environment variable overrides, which is useful when a collaborator does not
want to edit a local config file.

## Active PheWAS Scripts

1. `15_positive_lung_cultures_prior_year_pollution_landscape.R`

   Builds the analytic landscape file: every positive pulmonary culture row in
   the CLIF database linked to prior-year ZCTA PM2.5 and NO2.

   Current definitions:
   - Pulmonary specimens: `fluid_category` in `respiratory_tract` or
     `respiratory_tract_lower`.
   - Positive culture: `method_category == "culture"` and
     `organism_group`/`organism_category` is not `no_growth`.
   - Exposure year: `admission_year - 1`.
   - PM2.5: monthly ZCTA values annualized within ZIP/year.
   - NO2: annual ZCTA values.

   Main outputs:
   - `positive_lung_cultures_prior_year_pollution_<site>_<stamp>.csv`
   - `positive_lung_cultures_prior_year_pollution_coverage_<site>_<stamp>.csv`
   - `positive_lung_cultures_prior_year_pollution_organism_summary_<site>_<stamp>.csv`
   - `positive_lung_cultures_prior_year_pollution_year_summary_<site>_<stamp>.csv`

2. `17_plot_positive_lung_organism_phewas.R`

   Fits organism-specific logistic regression models and creates the main
   PheWAS-style plots. Each point is an `organism_category`, arranged by its
   dominant `organism_group`.

   Model:
   ```r
   organism_present ~ exposure_scaled + age_at_admission + sex_category +
     race_category + ethnicity_category + admission_year
   ```

   Current exposure scaling:
   - PM2.5: per 5 ug/m3.
   - NO2: per 10 ppb.

   Default minimum detections:
   - `MIN_ORGANISM_DETECTIONS=10`

   Main outputs:
   - `positive_lung_culture_organism_prior_year_pollution_models_<stamp>.csv`
   - `positive_lung_culture_organism_phewas_pm25_prior_year_<stamp>.png`
   - `positive_lung_culture_organism_phewas_no2_prior_year_<stamp>.png`

3. `16_plot_positive_lung_group_phewas.R`

   Optional companion screen at the `organism_group` level. This uses the same
   cohort, covariates, exposure scaling, and plot style, but the outcome is
   organism group present versus absent.

   Default minimum detections:
   - `MIN_GROUP_DETECTIONS=10`

## Run Order

From the repository root:

```bash
Rscript code/15_positive_lung_cultures_prior_year_pollution_landscape.R
Rscript code/17_plot_positive_lung_organism_phewas.R
```

Optional group-level companion:

```bash
Rscript code/16_plot_positive_lung_group_phewas.R
```

Buddy-test example without editing `config/config.json`:

```bash
CLIF_SITE_NAME=YOUR_SITE \
CLIF_TABLES_PATH=/path/to/CLIF/2.1.0 \
CLIF_FILE_TYPE=parquet \
ZCTA_EXPOSURE_DIR=data/exposome_zcta \
Rscript code/15_positive_lung_cultures_prior_year_pollution_landscape.R

MIN_ORGANISM_DETECTIONS=10 \
Rscript code/17_plot_positive_lung_organism_phewas.R
```

If a collaborator keeps the ZCTA exposure release somewhere else, set
`ZCTA_EXPOSURE_DIR` to that directory. It must contain:

1. `air_pollution_zcta_pm25_monthly_2005_2023.parquet`
2. `air_pollution_zcta_no2_annual_2005_2025.parquet`

## Configuration

Copy `config/config_template.json` to `config/config.json` and set:

1. `site_name`: short site label for output filenames.
2. `tables_path`: directory containing CLIF tables.
3. `file_type`: `parquet`, `csv`, `fst`, or `auto`.
4. `zcta_exposure_dir`: directory containing the ZCTA PM2.5 and NO2 parquet
   release files.

The table reader locates CLIF tables recursively under `tables_path` and accepts
filenames with or without the `clif_` prefix, as long as the base table name is
unique.

## Legacy Scripts

Older workflows live in subfolders so they do not distract from the current
PheWAS:

1. `legacy_county_level/`: original county-level aggregate pollution-microbe
   analyses.
2. `legacy_shrf_zcta/`: severe hypoxemic respiratory failure and ED-to-ICU
   early respiratory culture analyses.

Do not run these folders for buddy testing the current PheWAS unless the goal is
to reproduce those older analyses.

## Development Notes

The current code intentionally keeps exploratory outputs unsuppressed for local
signal finding. Before sharing outside an approved environment, review aggregate
release rules and add suppression if needed.
