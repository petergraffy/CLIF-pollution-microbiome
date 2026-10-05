# CLIF Pollution-Microbiome

This repository implements an organism-wide association study of acute residential air pollution and respiratory culture-detected organisms in CLIF 2.1. It uses a time-stratified case-crossover design and exports aggregate estimates for federated meta-analysis. Clinical cultures do not measure the complete sequencing-based respiratory microbiome.

## Quick start for sites

Use **R 4.4.2**, the version pinned in `renv.lock`. Run all commands from the repository root.

### 1. Clone the repository

```sh
git clone https://github.com/petergraffy/CLIF-pollution-microbiome.git
cd CLIF-pollution-microbiome
```

The clone includes the nationwide pollution, weather and ACS files (about 1.7 GB). Normal site runs require only R.

### 2. Restore the R environment

```sh
Rscript -e 'renv::restore(prompt = FALSE)'
```

The project environment activates automatically through `.Rprofile`. Restore once before the first run and again after an update changes `renv.lock`. The analysis pipeline does not install packages.

### 3. Configure your site

```sh
cp config/config_template.json config/config.json
```

Fill in `site_name`, `tables_path` and `file_type` (`parquet` for the current site pipeline). The local config is ignored by Git.

`site_timezone` defaults to `America/Chicago` (Central time). If your hospital is elsewhere, use its IANA time zone: `America/New_York` (Eastern), `America/Denver` (Mountain), `America/Los_Angeles` (Pacific), or `America/Phoenix` (Arizona). Use the hospital/data timestamp time zone, not your computer’s location; do not use abbreviations such as CST or EST. These names account for applicable daylight saving time.

All other settings use the shared protocol defaults, including bundled public-data paths and the first pilot’s disabled optional SOFA derivation.

### 4. Run the pipeline

```sh
Rscript code/00_run_pipeline.R
```

The pipeline checks the pinned packages, CLIF schemas and bundled data, creates a unique run, then performs cohort preparation, exposure linkage, Table 1/annual summaries, all models, diagnostics, local pooling and reporting. No Python setup or environmental download is needed.

### 5. Review and return results

The final message identifies the completed **`output/runs/<run_id>/`** folder. Open its `federated/report.html`, review coverage, culture mapping, missingness, sparse fits and calibration diagnostics, and check the empty `privacy_audit.csv`. After institutional disclosure review, return that one aggregate-only folder.

Patient-level working files remain in ignored `output/mwas/<run_id>/private/`. Do not return the working run, config, caches or logs. Small aggregate cells are not automatically suppressed; the export audit checks structured identifier fields and does not replace institutional review.

Optional synthetic check: `Rscript code/35_buddy_smoke_test.R`. To resume a failed run, use `Rscript code/00_run_pipeline.R --resume --run-id EXISTING_RUN_ID`. See [the buddy guide](docs/buddy_testing.md) and [script inventory](code/README.md) for details.

## Scientific definitions

- Hospital admission followed by first ICU entry within 24 elapsed hours; no age restriction.
- Primary respiratory culture window: first 48 ICU hours; 24/72-hour sensitivities.
- Source uses CLIF `fluid_category`, with culture methods only. Four respiratory categories are included; source-specific sensitivities are separate.
- Each named organism is an outcome, once per admission. All named taxa are attempted; no arbitrary 50/100-event cutoff.
- Daily PM2.5 per 5 µg/m³ and ozone per 10 ppb. Seven-day preceding mean is primary, with 3/14/28-day sensitivities. NO2 is excluded.
- Same-weekday/month/year reference dates, seven-day temperature/humidity splines and federal holidays; patient-cluster inference with conditional count-model checks.
- ACS neighborhood SES interactions, diagnosis-defined outcomes, age/sex/recorded-race and Charlson/modified-SOFA interactions.
- Table 1, annual site characteristics, culture practices, antibiotic timing and descriptive outcomes.
- Site estimates and compatible aggregate meta-analysis, with prespecified hypothesis families and failed-model tracking.

Severity and SES enter as pollution interactions. Their admission-level main effects cancel in self-matched models. SOFA totals require all six observed domains; optional physiology derivation is disabled in the default first buddy pilot. Enabling optional SOFA derivation is a coordinated protocol choice, not a required site configuration step.

## Bundled public inputs

[Public data and provenance](data/public/README.md) are included in the repository:

- Nationwide daily PM2.5, ozone, temperature and relative humidity, December 2017–December 2024, in monthly Parquet shards.
- Nationwide 2013–2017 ACS ZCTA SES indicators, the required source estimates/MOEs and Census variable metadata.
- SHA256 manifests and source release URLs.

Normal analysis reads these local files and makes no environmental downloads. National coverage supports another site without publishing a patient-derived ZIP subset. Only required daily columns are retained, without value rounding or imputation. Monthly files remain below GitHub's per-file size limit.

## Interpretation and sharing

This is a development protocol. Sparse results, culture selection, antibiotic exposure and post-entry severity require careful interpretation. Simulation has identified excess null rejection in some scenarios; technical success does not establish inferential validity.

Clinical records, local configuration, ZIP linkage lists, filtered site caches and generated outputs remain ignored and local. Site-derived aggregate outputs still require institutional disclosure review. Public inputs contain no patient information.
