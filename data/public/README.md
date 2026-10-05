# Bundled national public MWAS inputs

These files are public nationwide inputs, independent of CLIF patient records or residential ZIP selections. Ordinary site runs read them locally and do not download environmental data. Site-filtered historical caches remain ignored under `data/mwas_cache/` and must not be shared.

## Daily exposures

`exposures/{pm25,o3,weather}/` contains monthly Parquet shards for December 2017–December 2024. This supports every same-weekday/month/year candidate date in the 2018–2024 admission period with a maximum 28-day pre-date exposure window. No admission-day exposure is used.

Required columns:

| Product | Columns |
|---|---|
| PM2.5 | `zip`, `date`, `pm25_ug_m3`, `value_source`, `fill_distance_m` |
| Ozone | `zip`, `date`, `o3_ppb`, `value_source`, `fill_distance_m` |
| Weather | `zip`, `date`, `tmean_c`, `rhmean_pct` |

Source: [environment_transplant_survival releases](https://github.com/petergraffy/environment_transplant_survival/releases), tags `lghap-pm25-zcta-daily-v1`, `o3-zcta-daily-v1`, and `gridmet-zcta-daily-v1`. Every original annual asset is SHA256 verified against GitHub release metadata before conversion. Only unused columns are dropped; values retain source precision. ZIP strings are padded to five digits and dates stored as date32. Monthly files are sorted by ZIP/date and losslessly ZSTD compressed. No rounding, imputation, spatial selection or clinical filtering is performed. Duplicate ZIP/date records are rejected.

`exposures/manifest.json` records each shard's checksum, bytes, rows and unique ZCTA count, plus original source URLs/digests and transformations. Each product has annual `_source.json` provenance. Nation-wide scope follows each public source's available geography; it does not imply all residential ZIPs have coverage.

The frozen bundle has 255 monthly shards (about 1.70 GB total), each below 8.64 MB. Monthly shards keep individual files below 95 MB, avoiding GitHub's 100 MB file limit without Git LFS. The complete national bundle is substantially larger than a site-filtered cache. Git attributes preserve exact public-data bytes and LF source-code line endings across platforms. The runner validates checksums first; Arrow filters national files by required site ZIPs during local reading.

## ACS

`acs/2017/` contains all 32,989 published ZCTAs from the 2013–2017 ACS five-year summary:

- `zcta_ses.csv`: derived poverty, no-high-school education, unemployment, median household income and crowding indicators, with denominators.
- `acs_variables.csv`: required source estimates (`E`), margins of error (`M`) and available annotations. Summary-file bounded income values receive an explicit annotation.
- Five `B*_metadata.json` files: official Census variable descriptions.
- `manifest.json`: definitions, source URLs/checksums, bundled-file checksums, missing-value rules and uncertainty limitations.

Definitions: poverty `B17001_002/001`; education below high school among age 25+ `sum(B15003_002..016)/001`; unemployment `B23025_005/003`; median household income `B19013_001`; crowding `sum(B25014_005,006,007,011,012,013)/001`. Percentages use valid denominators. Negative sentinels, annotated medians, invalid ratios and zero denominators remain missing. Sampling uncertainty is retained but not propagated into model coefficients.

Source: [2017 ACS five-year Census summary files](https://www2.census.gov/programs-surveys/acs/summary_file/2017/). No Census API key is needed for the default summary-file rebuild.

## Offline validation and explicit maintainer rebuild

```bash
python3 code/22_cache_mwas_exposures.py
python3 code/26_cache_acs_zcta_ses.py
```

These commands only validate finished bundled files. Missing/corrupt files cause an error; the site runner never silently replaces them with network downloads.

To intentionally rebuild public data with the pinned Python environment:

```bash
.venv-buddy/bin/python code/22_cache_mwas_exposures.py --download
.venv-buddy/bin/python code/26_cache_acs_zcta_ses.py --download
```

The ACS maintainer raw-download cache is ignored at `data/mwas_cache/acs/`. Exposure temporary annual downloads are removed after sharding. Study periods outside 2018–2024 require a deliberate new data bundle and compatible protocol; extending the analysis dates alone is insufficient.

A full repository checkout includes the national inputs. `python3 code/36_package_buddy_source.py --include-public-data` creates an allowlisted source-and-public-data archive. Its default source-only archive requires a separate matching copy of this public bundle.
