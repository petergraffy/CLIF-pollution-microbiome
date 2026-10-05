# Local site configuration

Copy `config_template.json` to `config.json`, which is ignored by Git. Only three fields need site input:

| Field | What to enter |
|---|---|
| `site_name` | A unique site label, such as `UCMC` |
| `tables_path` | The local folder containing your CLIF 2.1 tables |
| `file_type` | `parquet` for the current site pipeline; use the extension without a leading dot |

`site_timezone` defaults to `America/Chicago` (Central time). If your hospital is elsewhere, use its IANA time zone: `America/New_York` (Eastern), `America/Denver` (Mountain), `America/Los_Angeles` (Pacific), or `America/Phoenix` (Arizona). Use the hospital/data timestamp time zone, not your computer’s location; do not use abbreviations such as CST or EST. These names account for applicable daylight saving time.

Actual calendar dates and correctly interpreted timestamp offsets are required.

Study dates, exposure windows, models, bundled public-data paths and ACS variables use the shared protocol defaults. Sites do not need to add these settings to their configuration. Standard first-pilot runs leave optional SOFA derivation disabled; Charlson remains part of the workflow. Enabling SOFA is a coordinated protocol choice requiring additional physiology tables, not a routine site setup step. Missing/skipped SOFA totals remain missing.

The six core tables are patient, hospitalization, adt, microbiology_culture, medication_admin_intermittent and hospital_diagnosis. Run `Rscript code/33_site_preflight.R` to check required fields. Ethnicity and admission type are optional descriptive fields.

Ordinary site runs use bundled national public inputs, make no environmental network requests and require no Census key. `CLIF_CONFIG_PATH` can select an alternate local configuration file. Never commit credentials, completed configurations, source tables or clinical linkage files.
