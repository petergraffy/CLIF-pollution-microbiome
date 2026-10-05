## Configuration

Copy `config_template.json` to `config.json` and update it for the local environment.

Required fields:

1. `site_name`: short site label used in output filenames.
2. `repo`: absolute path to this repository.
3. `tables_path`: absolute path to the CLIF table directory.
4. `file_type`: CLIF table file type, usually `parquet`, `csv`, or `fst`.
5. `zcta_exposure_dir`: path to the ZCTA air-pollution parquet release used by the current PheWAS workflow.

Example:

```json
{
    "site_name": "YOUR_SITE",
    "repo": "/path/to/CLIF-pollution-microbiome",
    "tables_path": "/path/to/CLIF/2.1.0",
    "file_type": "parquet",
    "zcta_exposure_dir": "/path/to/air_pollution_zcta_parquet"
}
```

The current PheWAS scripts expect these files in `zcta_exposure_dir`:

1. `air_pollution_zcta_pm25_monthly_2005_2023.parquet`
2. `air_pollution_zcta_no2_annual_2005_2025.parquet`

The committed `data/exposome_zcta` directory already contains these files and
is the default if `zcta_exposure_dir` is omitted.

The code locates CLIF tables recursively under `tables_path` and accepts filenames with or without the `clif_` prefix, as long as the base table name is unique. For example, `clif_hospitalization.parquet` and `hospitalization.parquet` are both valid.

Common environment variable overrides:

1. `CLIF_CONFIG_PATH`: alternate config JSON path.
2. `CLIF_SITE_NAME`: override `site_name`.
3. `CLIF_TABLES_PATH`: override `tables_path`.
4. `CLIF_FILE_TYPE`: override `file_type`; use `auto` to scan `csv`, `parquet`, and `fst`.
5. `ZCTA_EXPOSURE_DIR`: override `zcta_exposure_dir`.

Optional legacy field:

- `exposome_path`: county-year exposure directory used only by the older aggregate county-level scripts. It is not needed for the current PheWAS workflow.

The `.gitignore` file in this directory prevents `config.json` from being pushed to GitHub. Keep site-specific paths and credentials local.

For the acute MWAS, also set `site_timezone` to the hospital's IANA timezone
(e.g. `America/Chicago` for UCMC). Elapsed inclusion windows use timestamp
instants; daily exposure matching uses local admission dates. `tables_path`
should point to the specific CLIF release folder, not a parent containing
additional disease-specific exports. Daily exposures are cached separately
under `data/mwas_cache/`; see the acute workflow in `code/README.md`.

The federated SES extension uses `MWAS_ACS_DIR` (default
`data/mwas_cache/acs/2017`) for the national ACS cache. All sites must use the
same vintage and indicator definitions. The default public summary-file route
needs no key. For the optional API route, set `CENSUS_API_KEY` in the local
environment; never put it in committed configuration or share it in outputs.

The v4 buddy MWAS additionally requires the patient table (`patient_id`, `sex_category`, `race_category`) and `hospitalization.discharge_category`. Ethnicity and admission type are optional descriptive fields. The run writes private demographic/outcome linkage and aggregate Table 1/annual exports; see the [buddy guide](../docs/buddy_testing.md) for denominators and disclosure review.
