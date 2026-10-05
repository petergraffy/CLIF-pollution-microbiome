# Acute respiratory culture MWAS: buddy testing

This is a development protocol for a second CLIF 2.1 site. It tests portability, culture capture, inference and aggregate pooling. Clinical cultures measure detected organisms, not sequencing-based microbiome composition. No patient records leave the site.

## Site quick start

Use **R 4.4.2** and the repository root as the working directory. The ordinary site pipeline is R-only. Python/PyArrow are maintainer tools for rebuilding public data or packaging archives, not site-run prerequisites.

```sh
git clone https://github.com/petergraffy/CLIF-pollution-microbiome.git
cd CLIF-pollution-microbiome
Rscript -e 'renv::restore(prompt = FALSE)'
cp config/config_template.json config/config.json
```

Edit `site_name`, `tables_path` and `site_timezone`. Keep the bundled public-data paths and Parquet file type. `CLIF_CONFIG_PATH` can select an alternate local config. R packages are pinned in the single `renv.lock`; `.Rprofile` activates the environment automatically. Package restoration is a setup step and never happens inside the analysis pipeline. After changing the lockfile, restore again.

Required tables: patient, hospitalization, adt, microbiology_culture, medication_admin_intermittent and hospital_diagnosis. Confirm actual calendar dates, correct timestamp timezone interpretation and residential ZIP. The first pilot defaults to `derive_sofa=false`. Set it to `true` if labs, vitals, respiratory_support, patient_assessments and medication_admin_continuous are available. Complete six-domain SOFA totals are required; missing/skipped scores do not become normal scores. Six-hour SOFA is the primary acute modifier and 24-hour SOFA a sensitivity, both measured after ICU entry.

Optional synthetic R check before clinical data:

```sh
Rscript code/35_buddy_smoke_test.R
```

Run everything with:

```sh
Rscript code/00_run_pipeline.R
```

The runner creates a unique site/timestamp run ID, checks dependencies/schemas, validates national input SHA256 checksums, then runs the cohort, exposure matching, characteristics, models, count checks, local pooling, calibration, sparse audits, report and curated export. Study dates stay fixed at 2018–2024. Bundled daily inputs include December 2017 for the maximum lookback. Normal runs make no environmental network requests.

The final message prints **`output/runs/<run_id>/`**, the aggregate-only return folder. Open `federated/report.html`, review the diagnostics and empty `privacy_audit.csv`, and apply institutional small-cell/disclosure rules before returning that folder. `export_manifest.json` records file checksums and release-review status. No transmission happens automatically.

Patient-level working files and local diagnostics remain under ignored `output/mwas/<run_id>/`; never return that entire directory. The export uses an explicit allowlist, rejects identifier columns/JSON keys and excludes private files, raw linkage, local configuration and clinical logs. This structural audit is not a comprehensive de-identification guarantee.

After a failed stage:

```sh
Rscript code/00_run_pipeline.R --resume --run-id EXISTING_RUN_ID
```

Resume reuses the prepared cohort and recomputes analysis stages. Start a fresh run after clinical inputs or preparation definitions change. To prepare exports from an already completed working run without rerunning models:

```sh
MWAS_RUN_DIR=output/mwas/EXISTING_RUN_ID Rscript code/07_prepare_site_exports.R
```

An existing export destination is never overwritten. `MWAS_EXPORT_DIR` can select a new destination. `MWAS_EXPOSURE_CACHE` and `MWAS_ACS_DIR` remain advanced overrides for deliberate alternate public bundles. See [public-data provenance](../data/public/README.md).

## Frozen scientific definitions

- Eligibility: a hospital admission whose first ICU entry occurs 0–24 elapsed hours after hospital admission; no age restriction. Invalid residence ZIPs are excluded from exposure linkage and counted in cohort flow. The event is local hospital admission date, not collection date or ICU date.
- Specimens: `method_category == "culture"` and CLIF `fluid_category` in respiratory_tract, respiratory_tract_lower, nasopharynx_upperairway, oropharynx_tongue_oralcavity. `fluid_name` is audited, not used to redefine eligibility. The primary window is 0–48 elapsed ICU hours; 24 and 72 hours are sensitivities.
- Outcome: any named organism detected in an eligible culture, one event per organism per admission. All named taxa within 72 hours are registered and attempted; no 50/100-count filter. No-growth and nonspecific flora do not become named-taxon hypotheses. Sparse/failed models remain in the attempted registry. Organism taxonomy must be mapped consistently to CLIF categories across sites.
- Exposures: daily PM2.5 per 5 micrograms/m3 and ozone per 10 ppb. Complete means on days 1–3, 1–7 (primary), 1–14, 1–28 before admission. No admission-day exposure. NO2 is excluded. Weekday/month/year referents use the same residential ZIP proxy. Candidate rows and eligible cases are identical across all four windows within each pollutant, with complete weather and nonzero exposure variation in every window.
- Adjustment: natural splines (3 df each) of lag-1–7 temperature and relative humidity, plus US federal holiday indicator, fixed across durations. Conditional logistic models use patient-cluster covariance. Time-invariant patient severity or neighborhood SES cannot be fitted as a main effect within a case-only stratum. No inpatient-referent exclusion is applied in this primary protocol.
- ACS: poverty is the primary modifier, per 10 percentage points, centered at 20%. Education below high school (25+), unemployment (civilian labor force), income (log2 relative to $50,000), and crowding (>1 occupant/room) are secondary. Interactions estimate modification of the pollution OR; ACS uncertainty is retained but not propagated.
- Admission diagnosis: primary diagnosis flagged POA. Pneumonia/aspiration, obstructive airway, other respiratory and nonrespiratory outcomes use the existing ICD9/10 mapping. Missing/ambiguous/unmapped codes remain explicit. This is an outcome subgroup analysis, not proof of baseline ascertainment.

## Added checks (items 1–4)

1. **Shared exposure inference.** Conditional Poisson likelihood on identical ZIP/candidate-date groups reproduces conditional logistic point estimates. Separate uncertainty checks use Pearson quasi-dispersion (floored at 1) and calendar HAC with Bartlett weights through 28 calendar days. Calendar scores are projected off conditioned group intercepts and aggregated across all ZIPs on each date. This accounts for some shared temporal dependence; it does not guarantee valid inference for sparse taxa. SES interactions retain the original patient-cluster method and are not covered by this first count-model check. See the [conditional Poisson formulation](https://doi.org/10.1186/1471-2288-14-122).
2. **Selection companions.** All early ICU admissions, any respiratory culture within 48 ICU hours, and any named organism within 48 hours are modeled with the same exposure windows. These describe admission/testing patterns; they do not identify a causal organism-specific effect or eliminate conditioning-on-hospitalization bias.
3. **Clinical sensitivities.** At the seven-day exposure window, each taxon is tested for cultures within 24/72 hours, pulmonary categories only, upper-airway categories only, and detection before any documented inpatient antibacterial administration. The last sensitivity requires known medication history and excludes equal-time/unknown histories. Outpatient and referring-hospital antibiotic exposure remain unmeasured. Pulmonary `respiratory_tract` is only as specific as the site's category mapping.
4. **Simulation calibration.** Default 200 replicates per scenario, with 12/100/600 events, 7/28-day exposures, both pollutants, true OR 1 and 1.5, and independent versus site-wide AR(1) calendar shocks. Actual local exposure/weather series and candidate sets are used. Outputs include failed-fit rates, null false-positive rates, conditional coverage, OR 1.5 power and binomial Monte Carlo intervals. Failed fits count as nonrejections. This is conditional single-hypothesis calibration, not BH discovery power or a population selection simulation. Shared shocks marginalize a conditional effect, so bias need not be zero in that scenario. Use `MWAS_SIM_REPS=1000` for more precise calibration after the pilot; changing simulation settings does not change the fitted MWAS protocol.

## Demographic and severity effect modification

Protocol v4 adds seven-day pollution interactions fitted one modifier at a time. Age is continuous per 10 years, centered at 60; Charlson per point, centered at 2; complete six-hour modified SOFA per 2 points, centered at 6. Complete 24-hour SOFA uses the same scale as a sensitivity. Continuous interactions assume a linear change in the pollution log-OR over the modifier. Missing or invalid values are excluded; partial SOFA sums never substitute for complete totals.

Recorded sex compares Male with Female. Recorded race uses predeclared pairwise contrasts: Black or African American, Asian, American Indian or Alaska Native, Native Hawaiian or Other Pacific Islander, and Other, each compared separately with White. Each race model contains its category and the reference category; other categories and Unknown/Missing/Unmapped are excluded from that contrast. No rare categories are combined and every attempt is retained. These pairwise models have different category-specific analysis populations; do not treat their reference pollution estimates as a single shared estimate. There is no omnibus race test and the contrasts are not independent because they share a reference group. Race is a recorded social classification and should not be interpreted as a biological mechanism. Categories follow the [CLIF 2.1 dictionary](https://clif-icu.com/data-dictionary/data-dictionary-2.1.0).

Because SOFA is measured after ICU entry and can reflect exposure consequences and treatment, interpret its interactions as associations across events of different observed severity, rather than evidence that severity causally modifies the pollution effect.

The interaction coefficient is a log ratio of pollution ORs; its exponent is a ratio of pollution ORs. Main pollution terms are effects at the stated reference value/category and are not separate interaction hypotheses. Reference categories are parameterization choices, not claims of a normative group. Interactions are currently evaluated using patient-cluster covariance only; the baseline count checks and simulations do not validate these interactions. Source/model/ACS package protocol compatibility is advanced to v4, so v3 and v4 exports cannot be pooled together.

## Table 1 and annual site exports

`code/37_site_characteristics.R` writes `table1.csv` for reading and `table1_long.csv` for machine use. The primary cohort is early ICU admission within 24 hospital hours, valid ZIP, and any eligible respiratory culture within 48 ICU hours, including cultures without a named organism. Context columns show all early ICU admissions with valid ZIP, the primary cohort with matched exposure sets for both pollutants, and the subset with a named organism. These are overlapping populations, not independent comparison groups; no Table 1 p-values are computed. All statistics are admission-level, including repeat admissions; unique patients are shown separately. Continuous values show median/IQR and observed/total counts; categorical percentages use all admissions in the column, with unknown/missing shown explicitly.

Annual outputs contain site, local hospital admission year, and cohort-specific denominators for every study year, including zero-admission years:

- `site_year_characteristics.csv`: recorded CLIF hospitalizations, any ICU, all early ICU, early ICU with valid ZIP, primary cultures, matched primary cohort and named-organism subset; admissions/patients, culture coverage and episode rate, named detection yield, age/severity summaries, observed/missing outcome counts, hospital/ICU duration, death, hospice and documented IMV.
- `site_year_categories.csv`: sex, race, ethnicity, admission type, discharge disposition, and first-culture antibiotic status.
- `site_year_specimen_sources.csv`: primary-cohort culture episodes and admissions by CLIF fluid category, with episode percentages.
- `site_year_antibiotic_practices.csv`: primary culture episode antibiotic status, episode counts/percentages and admission counts. An admission can contribute episodes in more than one category.
- `site_year_organisms.csv`: named organism episode/admission counts in the primary cohort, assigned to admission year. Small cells require disclosure review.
- `characteristics_manifest.json`: denominator, timing, missingness and outcome definitions; `annual_site_trends.png` is a local descriptive plot generated with the report.

The hospital admission denominator is the available CLIF extract, not necessarily all admissions to the institution. Outcomes are assigned to admission year even when discharge occurs later. Death is `Expired`; hospice is separate. Mortality rates exclude incomplete admissions and missing/unmapped dispositions, with denominator counts retained. ICU duration is the union of recorded ICU intervals clipped to the hospital stay; invalid/incomplete intervals remain missing. IMV is any documented `imv` device during the admission; missing respiratory-support records remain unknown. Charlson is available only for early ICU admissions with valid ZIP, and SOFA derivation is restricted to those cultured within 72 hours. Descriptive missingness reflects these derivation scopes as well as clinical recording. Primary-cohort annual rows are the appropriate reference for culture-related severity summaries.

BH families are primary overall seven-day, primary poverty interaction, secondary SES interactions, diagnosis-defined outcomes, clinical sensitivities, selection companions, demographic modification, clinical modification (Charlson and six-hour SOFA), and 24-hour SOFA modification sensitivity. Nonprimary durations have separate sensitivity families, corrected jointly across 3/14/28-day windows. Corrections are separate by inference method and retain every attempted hypothesis, including failures. Do not select the smallest p-value across windows, methods or specimen definitions.

## Buddy review and aggregate return

First inspect cohort flow, exposure coverage, source capture and antibiotic status locally. Verify that respiratory categories reflect intended specimens at the site, medication vocabulary coverage is adequate, diagnosis POA capture is plausible, and culture windows use collection time. Reconcile count-model and conditional logistic estimates for successful fits. Review sparse estimates and simulation error rates before scientific interpretation. Technical success is not statistical validation. Development calibration has identified excess null rejection in some scenarios, especially for calendar HAC, and frequent failures with very sparse outcomes. The report flags null cells whose Monte Carlo interval lies above 5%; treat these methods as development checks until calibration is resolved.

Only release institutionally approved aggregate outputs. Suggested coordinator inputs are `federated/site_estimates.csv`, `federated/shared_exposure_estimates.csv` and `federated/protocol.json`. `simulation_calibration.csv`, Table 1 and annual site CSVs are useful diagnostics after local review. Reference/contrast definitions and interaction covariance are preserved for central pooling; `modifier_registry.csv` and `modifier_capture.csv` document model scales and missingness. `modifier_case_support.csv` gives admission- and patient-level category/score support for every attempted interaction; `modifier_candidate_diagnostics.csv` audits local FDR-passing interactions with model-based uncertainty and leave-one-patient-out refits. A large total event count does not guarantee adequate category-specific support. Counts and organism labels can still require institutional small-cell suppression. Apply the site's approved disclosure rules; preserve every attempted hypothesis row with an explicit suppressed status and remove its coefficient/SE as needed. Never send `private/`, source tables, ZIP linkage files, exposure lookup subsets, local config, clinical logs or the entire output/cache directory. The buddy archives contain none of these clinical or site-specific files. Public data can be included explicitly with `--include-public-data`.

The coordinator pools approved CSVs from each site together:

```bash
Rscript code/28_pool_federated_mwas.R output/mwas/pooled \
  /approved/site_A/site_estimates.csv /approved/site_A/shared_exposure_estimates.csv \
  /approved/site_B/site_estimates.csv /approved/site_B/shared_exposure_estimates.csv
```

Protocol hashes must agree (code, scientific definitions, common ACS data and model package versions), and duplicate site/model rows are rejected. Fixed-effect inverse-variance pooling is primary; Q, I-squared and random-effect summaries are descriptive with few sites. Patient counts sum across sites and do not establish globally unique patients. A second site's real run and data-mapping review are still required before calling the workflow multisite validated.

## Source bundle

Use a full repository checkout for a data-complete buddy run. `python3 code/36_package_buddy_source.py --include-public-data` produces a source-and-national-public-data archive in `dist/`. Without that flag the archive is source-only and requires the matching `data/public/` bundle from the repository. Both modes use allowlists and exclude patient data, local config, site caches, outputs, unrelated drafts and Git history. Share identical source/data versions to keep protocol hashes compatible. The v5 protocol also hashes the national exposure manifest and ACS indicator file, so different public bundles cannot silently pool. Older exploratory scripts were removed from `code/`; Git history preserves them.
