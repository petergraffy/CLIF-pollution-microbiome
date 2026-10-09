# MWAS results and local working files

Only this README is tracked in Git. All generated outputs are ignored and must be shared through approved aggregate exports, not repository commits.

The recommended pipeline prepares `output/runs/<run_id>/`, an aggregate-only return folder with `privacy_audit.csv`, `export_manifest.json` and `federated/report.html`. Review it under institutional disclosure rules, then share only that completed folder. Export destinations are never overwritten.


Current runs are written to ignored `output/mwas/<run ID>/` directories. Older exploratory outputs may remain locally but are not used by the current workflow.

- `private/`: admission-level cohort, culture episodes, organism detections, exposure matching, diagnosis and demographic/outcome linkage. Never distribute.
- Cohort, specimen, antibiotic and exposure QC: review locally for mapping and coverage.
- `federated/site_estimates.csv` and `shared_exposure_estimates.csv`: aggregate model terms, uncertainty, counts and attempted/failed status.
- `federated/protocol.json`: version, scientific definitions and compatibility hashes.
- `federated/table1.csv`, `table1_long.csv`: primary cohort and context populations, admission denominators and missingness.
- `federated/site_year_characteristics.csv`, `site_year_categories.csv`, `site_year_specimen_sources.csv`, `site_year_antibiotic_practices.csv`, `site_year_organisms.csv`: descriptive site trends with explicit denominators.
- Modifier registries/capture/support and candidate diagnostics: missingness, subgroup support and stability.
- Shared-exposure and simulation diagnostics: development checks; successful computation does not establish calibrated inference.
- Local pool, HTML report and figures: inspect before scientific interpretation or release.

Only institutionally approved aggregates may be shared. Small cells, unsuppressed organism counts, clinical logs, ZIP linkage files, private extracts and site-filtered caches must remain local unless separately authorized. Public nationwide pollution/weather/ACS inputs are kept separately in tracked `data/public/` and contain no clinical records.
