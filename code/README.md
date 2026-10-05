# Current acute MWAS scripts

Run from the repository root. The complete installation, scientific definitions and disclosure rules are in [the buddy guide](../docs/buddy_testing.md).

| Script | Purpose |
|---|---|
| `00_run_pipeline.R` | Recommended single-command site entry point |
| `07_prepare_site_exports.R` | Curated aggregate-only return folder, structural privacy audit and checksums |
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

Use R 4.4.2 and the repository root. Restore packages separately, configure your site, then run:

```sh
Rscript -e 'renv::restore(prompt = FALSE)'
cp config/config_template.json config/config.json
# Fill in site_name, tables_path and file_type.
# Change site_timezone only if the hospital is outside Central time.
Rscript code/00_run_pipeline.R
```

The runner validates bundled public files in R and prints the completed `output/runs/<run_id>/` aggregate-only folder. Review `federated/report.html`, `privacy_audit.csv` and small cells before sharing. Private working files remain under `output/mwas/<run_id>/private/`.

`Rscript code/35_buddy_smoke_test.R` tests the R workflow on synthetic inputs, including the main entry point and export. Python scripts are optional maintainer tools; rebuilding data needs the pinned PyArrow dependency, while normal site analysis requires no Python.

Public inputs and optional physiology derivation use shared defaults; sites do not need to configure these. See [config guidance](../config/README.md) for hospital time zones. To resume: `Rscript code/00_run_pipeline.R --resume --run-id EXISTING_RUN_ID`. Script 34 remains the underlying orchestrator; users need only the 00 entry point.

To export an already completed working run, set `MWAS_RUN_DIR` and run `Rscript code/07_prepare_site_exports.R`. Source/public-data archive tools are described in [the buddy guide](../docs/buddy_testing.md).

## Coordinator

```bash
Rscript code/28_pool_federated_mwas.R output/mwas/pooled \
  /path/to/approved/site1_estimates.csv /path/to/approved/site2_estimates.csv
```

Include approved shared-exposure estimate files as additional arguments to pool those methods. Pooling rejects protocol mismatches and duplicate site/model records. Patient-level records are never coordinator inputs.
