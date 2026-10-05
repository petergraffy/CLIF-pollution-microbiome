# Current acute MWAS scripts

Run from the repository root. The complete installation, scientific definitions and disclosure rules are in [the buddy guide](../docs/buddy_testing.md).

| Script | Purpose |
|---|---|
| `21_prepare_acute_mwas.R` | Early-ICU cohort, respiratory culture episodes, named detections, antibiotic history, Charlson and optional modified SOFA |
| `22_cache_mwas_exposures.py` | Offline checksum validation of bundled nationwide exposures; explicit `--download` maintainer rebuild |
| `23_run_acute_mwas.R` | Calendar reference dates, local public exposure linkage, exact 3/7/14/28-day means and coverage QC; prepares exposures only |
| `26_cache_acs_zcta_ses.py` | Offline validation of national ACS bundle; explicit `--download` derivation/rebuild |
| `27_run_federated_mwas.R` | All named-organism models, diagnosis/clinical/selection analyses, SES/demographic/severity interactions, aggregate exports and protocol hash |
| `28_pool_federated_mwas.R` | Compatible site pooling, covariance, coverage, hypothesis-family FDR and power diagnostics |
| `29_report_federated_mwas.R` | Local HTML report, duration plots and annual site trends |
| `30_audit_duration_candidates.R` | Sparse duration candidate uncertainty and leave-one-patient-out diagnostics |
| `31_shared_exposure_checks.R` | Conditional Poisson quasi-dispersion and calendar-HAC inference checks |
| `32_calibrate_mwas.R` | Known-null/OR 1.5 simulations using local exposure/reference-date patterns |
| `33_site_preflight.R` | Core CLIF schema, timezone and dependency validation |
| `34_run_buddy_site.R` | Full local workflow; explicit `--resume` |
| `35_buddy_smoke_test.R` | Unit checks, synthetic full pipeline and two-site pooling tests |
| `36_package_buddy_source.py` | Allowlisted source archive; `--include-public-data` adds only the national public bundle |
| `37_site_characteristics.R` | Private demographic/outcome linkage; aggregate Table 1 and annual characteristics |
| `38_audit_modifier_support.R` | Modifier-specific support and sparse interaction diagnostics |

Number gaps preserve familiar script names and references. Legacy county-level/SHRF/prior-year workflows, plotting-only grant scripts, templates and superseded summary/power scripts have been removed. Git history preserves earlier versions.

## Run

```bash
export MWAS_RUN_ID=MY_SITE_first_pass
export MWAS_PYTHON="$(pwd)/.venv-buddy/bin/python"
Rscript code/35_buddy_smoke_test.R
Rscript code/34_run_buddy_site.R
```

The runner validates `data/public/exposures/manifest.json` and `data/public/acs/2017/manifest.json` offline. Environmental files are already downloaded. Public-data rebuild commands and definitions are in [data/public/README.md](../data/public/README.md). Source reads use Arrow filters to load only required ZIPs into local memory; those selections remain private.

To enable six-/24-hour modified SOFA, set `MWAS_SKIP_SEVERITY=0` before a fresh run. `MWAS_EXPOSURE_CACHE` and `MWAS_ACS_DIR` deliberately override bundled paths. Generated output is `output/mwas/<run ID>/`; do not publish its private data or unreviewed aggregates.

## Coordinator

```bash
Rscript code/28_pool_federated_mwas.R output/mwas/pooled \
  /path/to/approved/site1_estimates.csv /path/to/approved/site2_estimates.csv
```

Include approved shared-exposure estimate files as additional arguments to pool those methods. Pooling rejects protocol mismatches and duplicate site/model records. Patient-level records are never coordinator inputs.
