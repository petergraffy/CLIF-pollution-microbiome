# CLIF Pollution-Microbiome

This repository implements an organism-wide association study of acute residential air pollution and respiratory culture-detected organisms in CLIF 2.1. It uses a time-stratified case-crossover design and exports aggregate estimates for federated meta-analysis. Clinical cultures do not measure the complete sequencing-based respiratory microbiome.

## Start here

1. Follow the [buddy-site installation and analysis guide](docs/buddy_testing.md).
2. Run `Rscript code/35_buddy_smoke_test.R` using synthetic data.
3. Copy `config/config_template.json` to the ignored `config/config.json` and set site name, CLIF table path and timezone.
4. Set a new `MWAS_RUN_ID`, then run `Rscript code/34_run_buddy_site.R`.

The [script inventory](code/README.md) lists only the current workflow. Older county-level, prior-year PheWAS, grant-figure and template scripts were removed from `code/`; they remain recoverable in Git history.

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

Severity and SES enter as pollution interactions. Their admission-level main effects cancel in self-matched models. SOFA totals require all six observed domains; optional physiology derivation is disabled in the default first buddy pilot. Set `MWAS_SKIP_SEVERITY=0` to enable it.

## Bundled public inputs

[Public data and provenance](data/public/README.md) are included in the repository:

- Nationwide daily PM2.5, ozone, temperature and relative humidity, December 2017–December 2024, in monthly Parquet shards.
- Nationwide 2013–2017 ACS ZCTA SES indicators, the required source estimates/MOEs and Census variable metadata.
- SHA256 manifests and source release URLs.

Normal analysis reads these local files and makes no environmental downloads. National coverage supports another site without publishing a patient-derived ZIP subset. Only required daily columns are retained, without value rounding or imputation. Monthly files remain below GitHub's per-file size limit.

## Interpretation and sharing

This is a development protocol. Sparse results, culture selection, antibiotic exposure and post-entry severity require careful interpretation. Simulation has identified excess null rejection in some scenarios; technical success does not establish inferential validity.

Clinical records, local configuration, ZIP linkage lists, filtered site caches and generated outputs remain ignored and local. Site-derived aggregate outputs still require institutional disclosure review. Public inputs contain no patient information.
