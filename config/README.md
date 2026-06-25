## Configuration

Copy `config_template.json` to `config.json` and update it for the local environment.

Required fields:

1. `site_name`: short site label used in output filenames.
2. `repo`: absolute path to this repository.
3. `tables_path`: absolute path to the CLIF table directory.
4. `file_type`: CLIF table file type, usually `parquet`, `csv`, or `fst`.
5. `exposome_path`: absolute path to county-year pollution exposure files used by the original county-level workflow.
6. `zcta_exposure_dir`: path to the ZCTA air-pollution parquet release used by the severe hypoxemic respiratory failure workflow.

Example:

```json
{
    "site_name": "YOUR_SITE",
    "repo": "/path/to/CLIF-pollution-microbiome",
    "tables_path": "/path/to/CLIF/2.1.0",
    "file_type": "parquet",
    "exposome_path": "/path/to/exposome",
    "zcta_exposure_dir": "/path/to/air_pollution_zcta_parquet"
}
```

The SHRF/ZCTA scripts expect these files in `zcta_exposure_dir`:

1. `air_pollution_zcta_pm25_monthly_2005_2023.parquet`
2. `air_pollution_zcta_o3_monthly_2005_2023.parquet`
3. `air_pollution_zcta_no2_annual_2005_2025.parquet`

The code locates CLIF tables recursively under `tables_path` and accepts filenames with or without the `clif_` prefix, as long as the base table name is unique. For example, `clif_hospitalization.parquet` and `hospitalization.parquet` are both valid.

Common environment variable overrides:

1. `CLIF_CONFIG_PATH`: alternate config JSON path.
2. `CLIF_SITE_NAME`: override `site_name`.
3. `CLIF_TABLES_PATH`: override `tables_path`.
4. `CLIF_FILE_TYPE`: override `file_type`; use `auto` to scan `csv`, `parquet`, and `fst`.
5. `ZCTA_EXPOSURE_DIR`: override `zcta_exposure_dir`.

The `.gitignore` file in this directory prevents `config.json` from being pushed to GitHub. Keep site-specific paths and credentials local.
