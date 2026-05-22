## Configuration

Copy `config_template.json` to `config.json` and update it for the local environment.

Required fields:

1. `site_name`: short site label used in output filenames.
2. `repo`: absolute path to this repository.
3. `tables_path`: absolute path to the CLIF table directory.
4. `file_type`: CLIF table file type, usually `parquet`, `csv`, or `fst`.
5. `exposome_path`: absolute path to county-year pollution exposure files.

Example:

```json
{
    "site_name": "UCMC",
    "repo": "/path/to/CLIF-pollution-microbiome",
    "tables_path": "/path/to/CLIF/2.1.0",
    "file_type": "parquet",
    "exposome_path": "/path/to/exposome"
}
```

The `.gitignore` file in this directory prevents `config.json` from being pushed to GitHub. Keep site-specific paths and credentials local.
