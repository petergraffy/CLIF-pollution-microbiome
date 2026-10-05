## Code Directory

This directory contains the executable R workflow for the current CLIF
pollution-microbiome project. The active workflow is now the positive lung
culture organism-wide pollution association screen, or PheWAS-style workflow.

Scripts read site-specific paths from `config/config.json` through
`utils/config.R` and `utils/clif_io.R`. They can also be buddy tested with
environment variable overrides, which is useful when a collaborator does not
want to edit a local config file.

## Active PheWAS Scripts

1. `15_positive_lung_cultures_prior_year_pollution_landscape.R`

   Builds the analytic landscape file: every positive pulmonary culture row in
   the CLIF database linked to prior-year ZCTA PM2.5 and NO2.

   Current definitions:
   - Pulmonary specimens: `fluid_category` in `respiratory_tract` or
     `respiratory_tract_lower`.
   - Positive culture: `method_category == "culture"` and
     `organism_group`/`organism_category` is not `no_growth`.
   - Exposure year: `admission_year - 1`.
   - PM2.5: monthly ZCTA values annualized within ZIP/year.
   - NO2: annual ZCTA values.

   Main outputs:
   - `positive_lung_cultures_prior_year_pollution_<site>_<stamp>.csv`
   - `positive_lung_cultures_prior_year_pollution_coverage_<site>_<stamp>.csv`
   - `positive_lung_cultures_prior_year_pollution_organism_summary_<site>_<stamp>.csv`
   - `positive_lung_cultures_prior_year_pollution_year_summary_<site>_<stamp>.csv`

2. `17_plot_positive_lung_organism_phewas.R`

   Fits organism-specific logistic regression models and creates the main
   PheWAS-style plots. Each point is an `organism_category`, arranged by its
   dominant `organism_group`.

   Model:
   ```r
   organism_present ~ exposure_scaled + age_at_admission + sex_category +
     race_category + ethnicity_category + admission_year
   ```

   Current exposure scaling:
   - PM2.5: per 5 ug/m3.
   - NO2: per 10 ppb.

   Default minimum detections:
   - `MIN_ORGANISM_DETECTIONS=10`

   Main outputs:
   - `positive_lung_culture_organism_prior_year_pollution_models_<stamp>.csv`
   - `positive_lung_culture_organism_phewas_site_summary_<stamp>.csv`
   - `positive_lung_culture_organism_phewas_modeled_organisms_<stamp>.csv`
   - `positive_lung_culture_organism_phewas_pm25_prior_year_<stamp>.png`
   - `positive_lung_culture_organism_phewas_no2_prior_year_<stamp>.png`

3. `16_plot_positive_lung_group_phewas.R`

   Optional companion screen at the `organism_group` level. This uses the same
   cohort, covariates, exposure scaling, and plot style, but the outcome is
   organism group present versus absent.

   Default minimum detections:
   - `MIN_GROUP_DETECTIONS=10`

   Main outputs:
   - `positive_lung_culture_group_prior_year_pollution_models_<stamp>.csv`
   - `positive_lung_culture_group_phewas_site_summary_<stamp>.csv`
   - `positive_lung_culture_group_phewas_modeled_groups_<stamp>.csv`

## Run Order

From the repository root:

```bash
Rscript code/15_positive_lung_cultures_prior_year_pollution_landscape.R
Rscript code/17_plot_positive_lung_organism_phewas.R
```

Optional group-level companion:

```bash
Rscript code/16_plot_positive_lung_group_phewas.R
```

Buddy-test example without editing `config/config.json`:

```bash
CLIF_SITE_NAME=YOUR_SITE \
CLIF_TABLES_PATH=/path/to/CLIF/2.1.0 \
CLIF_FILE_TYPE=parquet \
ZCTA_EXPOSURE_DIR=data/exposome_zcta \
Rscript code/15_positive_lung_cultures_prior_year_pollution_landscape.R

MIN_ORGANISM_DETECTIONS=10 \
Rscript code/17_plot_positive_lung_organism_phewas.R
```

If a collaborator keeps the ZCTA exposure release somewhere else, set
`ZCTA_EXPOSURE_DIR` to that directory. It must contain:

1. `air_pollution_zcta_pm25_monthly_2005_2023.parquet`
2. `air_pollution_zcta_no2_annual_2005_2025.parquet`

## Configuration

Copy `config/config_template.json` to `config/config.json` and set:

1. `site_name`: short site label for output filenames.
2. `tables_path`: directory containing CLIF tables.
3. `file_type`: `parquet`, `csv`, `fst`, or `auto`.
4. `zcta_exposure_dir`: directory containing the ZCTA PM2.5 and NO2 parquet
   release files.

The table reader locates CLIF tables recursively under `tables_path` and accepts
filenames with or without the `clif_` prefix, as long as the base table name is
unique.

## Legacy Scripts

Older workflows live in subfolders so they do not distract from the current
PheWAS:

1. `legacy_county_level/`: original county-level aggregate pollution-microbe
   analyses.
2. `legacy_shrf_zcta/`: severe hypoxemic respiratory failure and ED-to-ICU
   early respiratory culture analyses.

Do not run these folders for buddy testing the current PheWAS unless the goal is
to reproduce those older analyses.

## Development Notes

The current code intentionally keeps exploratory outputs unsuppressed for local
signal finding. Before sharing outside an approved environment, review aggregate
release rules and add suppression if needed.

## Acute PM2.5 / ozone MWAS (current development workflow)

This workflow studies organism-positive admissions among hospitalizations entering
an ICU within 24 elapsed hours of hospital admission. No additional age restriction
is applied. Respiratory culture windows are 0–24, 0–48 (primary), and 0–72 hours
from first ICU entry; pre-ICU specimens are inventoried but excluded. Exposure is
residential ZCTA pollution on days 1–7 before the local hospital admission date.
NO2 is excluded. Dates must be actual calendar dates.

Run from the repository root:

```bash
MWAS_RUN_ID=UCMC_first_pass Rscript --vanilla code/21_prepare_acute_mwas.R
python3 code/22_cache_mwas_exposures.py --run-dir output/mwas/UCMC_first_pass
MWAS_RUN_DIR=output/mwas/UCMC_first_pass Rscript --vanilla code/23_run_acute_mwas.R
Rscript --vanilla tests/test_mwas.R
```

Set `tables_path` and `site_timezone` in the ignored `config/config.json`.
Overrides include `CLIF_TABLES_PATH`, `MWAS_TIMEZONE`, `MWAS_START_DATE`,
`MWAS_END_DATE`, `MWAS_RUN_ID`, `MWAS_RUN_DIR`, `MWAS_EXPOSURE_CACHE`, and
`MWAS_MIN_EVENTS` (default 100 positive admissions per model). Default dates are
2018–2024. Exposure cache downloads default to 2017–2024; change the downloader's
`--first-year`/`--last-year` when changing the clinical date range. It needs Python
with PyArrow plus `curl`. The downloader fetches only public URLs; the local ZIP
list is never uploaded. Full yearly assets are checksummed and deleted after
retaining the needed ZCTAs. Release and subset metadata are saved with each run.

`output/mwas/` and `data/mwas_cache/` are ignored in Git. The `private/` run
subdirectory holds identifiers, dates, culture histories, and matched exposure
rows. Other run files contain unsuppressed LOCAL aggregate QC/results. None are
approved for external release. No patient data or site config is committed.

Models compare admission-day exposure histories with other weekdays of the same
month/year at the same residential ZCTA. Conditional logistic regression uses
Efron likelihood (one case per stratum), patient-cluster robust SE, natural
splines of the lag-1–7 mean temperature and relative humidity (3 df each), and
federal holidays. Main-effect Charlson, SOFA, and post-admission antibiotics are
not inserted as case-only covariates. FDR is corrected across named organisms and
both pollutants together within each analysis/window. The 48h all-cases family is
primary; 24h/72h, documented pre-antibacterial specimens, and excluding known
inpatient referents are secondary. Any-culture and any-named-organism companion
outcomes have a separate family. Failed/warning fits are retained as diagnostics
and do not produce reportable estimates.

Important first-pass definitions:

- CLIF lacks a specimen ID. Culture episode proxy: hospitalization, order and
  collection times, fluid category, method name, and LOINC. Serial result updates
  are consolidated per organism category. Simultaneous tests can be collapsed;
  culture revisions need further site review.
- Explicit negative text overrides category mapping (e.g. UCMC maps “no legionella
  isolated” to `legionella_sp`). Generic flora/bacteria/yeast outcomes are excluded
  from the named-organism screen. Taxon categories may still include genus-level
  results, not only species. Positive detection does not establish infection.
- Antibacterial administration uses the pinned CLIF mCIDE mapping, excludes
  antiviral/antifungal agents, and requires `mar_action_group=administered` before
  each specimen. Equal timestamps are uncertain. The primary flag includes all
  antibacterial routes and gut-directed agents; route-specific/spectrum-specific
  sensitivities remain future work. “None documented” does not establish absence
  of outpatient or transferring-hospital treatment. Observing any medication row
  establishes record availability, not completeness of antibiotic capture.
- Charlson uses POA ICD-9/10 Quan mappings and original Charlson weights with
  disease hierarchy. Missing diagnosis capture remains NA.
- Modified SOFA summaries use 0–6h and 0–24h after first ICU entry. These are early
  severity measures, not pre-admission baseline scores. Renal scoring is
  creatinine-only; respiratory scoring requires measured PaO2 paired with FiO2
  and support within one hour. Missing components remain missing; a full total
  requires all six domains. Partial sums are explicitly labeled. Sedation effects
  on GCS and unrecognized vasopressor units need additional validation.
- The current implementation is Stage 1. DLNM characterization, formal severity
  interactions, and external-site replication are not yet implemented.

A separately corrected exploratory 48h screen lowers the event threshold to 50.
These initial count cutoffs are provisional, not power-derived eligibility rules.

`Rscript code/25_audit_mwas_power_specimens.R` audits the existing 48h cohort
without changing its results. It writes original UCMC `fluid_name` labels and
design-based power scenarios for every named organism and both pollutants.
The current culture filter is `method_category == culture` and `fluid_category`
in `respiratory_tract`, `respiratory_tract_lower`, `nasopharynx_upperairway`, or
`oropharynx_tongue_oralcavity`. UCMC's original fluid names often describe the
assay rather than an anatomical specimen source; do not infer sputum/BAL from
these names. Resolve discordant labels before defining a lower-airway analysis.

Power uses null conditional likelihood information after within-set centering
and adjustment for the same weather spline bases and holiday term. Scenarios
cover OR 1.05, 1.1, 1.2, 1.5, and 2.0 per 5 ug/m3 PM2.5 or 10 ppb ozone, with 80%
target power at nominal alpha .05 and conservative Bonferroni alpha over all
named organisms and both pollutants. This is a local normal approximation,
assuming independent matched sets and zero nuisance coefficients; it is not
observed-effect post hoc power, nor BH discovery power. Required event counts
extrapolate each organism's existing exposure distribution. Validate by
simulation accounting for weather effects, repeated patients, and shared
calendar exposures before replacing eligibility rules. No final effect-size
target or power-based cutoff has yet been selected.

## Federated MWAS with ACS SES and admission diagnoses

The current extension attempts every named organism at 48h, retaining explicit
failed/insufficient-design rows. It does not select taxa or site estimates by
p-value. A model must have more independent patient clusters and matched events
than its within-set design rank plus one, identifiable exposure terms, no fit
warnings, and finite positive target covariance. These are computational
safeguards, **not** validation of sparse-sample inference or adequate power.
The original count-filtered outputs are retained as first-pass provenance.

Exposure-duration sensitivity now uses days **1–3, 1–7, 1–14, and 1–28** before
the local admission date. Seven days remains primary. Preparation caches daily
lags through 28 and computes a mean only with complete daily coverage; admission
day is excluded. Within each pollutant, site models use identical case and
referent rows across all four windows, with at least one referent and nonzero
variation in every window. SES and diagnosis restrictions are then applied to
that common population. PM2.5 and ozone populations can differ. Aggregate
`exposure_window_qc.csv` compares window-specific availability with the common
population. Temperature and humidity adjustment stay at lag 1–7 (3 df each),
and holidays and clinical definitions stay fixed: this isolates the pollution
window comparison rather than simultaneously changing weather adjustment.
Longer-window weather confounding remains a limitation of this sensitivity.

The 3-, 14-, and 28-day overall MWAS results form one joint
`exposure_duration_sensitivity` BH family across organisms and both pollutants.
Nonprimary diagnosis, poverty-interaction, and other-SES-interaction results
have their own corresponding duration sensitivity families, each combining all
three nonprimary windows. Seven-day families remain separate. Site exports,
pooling, interaction covariance, and power outputs include `exposure_window`;
the updated protocol prevents mixing the initial seven-day-only package with
this one. Window estimates are correlated, and comparing their p-values is not
a formal test of differing effects. The 28-day average can have much less
within-month exposure contrast.

Run the existing cohort/exposure stages once, then the new site workflow:

```bash
MWAS_RUN_ID=YOUR_SITE_run Rscript code/21_prepare_acute_mwas.R
python3 code/22_cache_mwas_exposures.py --run-dir output/mwas/YOUR_SITE_run
MWAS_RUN_DIR=output/mwas/YOUR_SITE_run MWAS_EXPOSURES_ONLY=1 Rscript code/23_run_acute_mwas.R
python3 code/26_cache_acs_zcta_ses.py
MWAS_RUN_DIR=output/mwas/YOUR_SITE_run Rscript code/27_run_federated_mwas.R
```

ACS defaults to national 2013–2017 five-year estimates, a fixed baseline preceding
the 2018–2024 admission period. Official public summary files are downloaded by
default, using the Census sequence lookup and national geography file. This does
not require an API key. `--source api` instead reads `CENSUS_API_KEY` from the
environment; credentials are not recorded in source manifests. Requests always
download national ZCTAs, never clinical ZIP lists. Estimates, margins of error,
metadata, source URLs and checksums are retained under the ignored
`data/mwas_cache/acs/2017/`. `MWAS_ACS_DIR` selects a local cache directory.
Every site must use the same vintage and national derived file.

| Indicator | ACS detailed table definition | Interaction scale / reference |
|---|---|---|
| Poverty (primary SES modifier) | B17001_002 / B17001_001 | 10 percentage points; reference 20% |
| No high school diploma, age 25+ | Sum B15003_002–016 / B15003_001 | 10 points; reference 20% |
| Unemployment, civilian labor force | B23025_005 / B23025_003 | 5 points; reference 5% |
| Median household income | B19013_001 | Doubling; reference $50,000 in 2017 dollars |
| Crowding, >1 occupants/room | Sum B25014_005–007 and 011–013 / B25014_001 | 5 points; reference 2% |

These are separate area-level indicators, not an invented deprivation composite
or individual SES measurements. No SES main effect is included: it is constant
within a matched admission and cancels from the conditional likelihood.
`pollution_ses` is the interaction coefficient; its exponent is a **ratio of
pollution odds ratios**, not the effect of poverty itself. `pollution` in an SES
model is the pollution effect at the stated reference. Main/interaction
covariance is exported and preserved in pooling for joint contrasts.

ZIP-to-ZCTA linkage uses exact five-digit codes as a proxy. It does not imply
postal ZIP and Census ZCTA boundaries are identical or reconcile 2010/2020 ZCTA
boundary changes. Missing geographic matches, negative sentinel values,
annotated/bounded income estimates, zero denominators, and invalid ratios remain
missing. No patient is excluded from the overall analysis for missing SES; SES
models use complete linked data. The fixed baseline may miss neighborhood
changes, and ACS sampling uncertainty is not yet propagated into interaction
estimates. Raw margins of error support later sensitivity work.

Diagnosis strata use `diagnosis_primary == 1` and `poa_present == 1`. Codes are
normalized for punctuation and classified using their ICD version:

| Group | ICD10 | ICD9 |
|---|---|---|
| Pneumonia/aspiration | J12–J18, J69 | 480–486, 507 |
| Obstructive airway | J41–J46 | 491, 492, 493, 496 |
| Other respiratory | Remaining J codes | Remaining 460–519 |
| Nonrespiratory | Other recognized ICD10 codes | Other recognized ICD9 codes |

These are explicit broad code groups, not a validated CCSR phenotyping package.
Influenza pneumonia is retained under other respiratory in this first definition.
Missing POA, absent primary diagnosis, unmapped codes, primary-not-POA, and
conflicting primary groups remain explicit QC categories, not nonrespiratory.
Primary discharge coding does not establish what clinicians knew on arrival;
sepsis may be principal despite a pulmonary presentation. These strata define
selected secondary outcomes, not adjustments for baseline confounding.

`federated/site_estimates.csv` uses an explicit aggregate schema: site/protocol,
organism/pollutant/analysis/term, coefficients, patient-cluster standard errors,
main/interaction covariance, event/patient/referent counts, design rank,
null information, and model status. It contains no patient IDs, dates, ZIPs,
clinical rows, or raw warning/error strings. Other files in `federated/` provide
aggregate linkage/diagnosis QC and protocol/provenance. Files stay local; apply
the site's aggregate-release policy before sharing. `private/` is never a site
deliverable. No upload or message to another site is performed by these scripts.

Central pooling uses site outputs only:

```bash
Rscript code/28_pool_federated_mwas.R output/mwas/pooled SITE_A/site_estimates.csv SITE_B/site_estimates.csv
```

The coordinator rejects differing protocol hashes (including code, study dates,
ACS vintage/data, specimen and model settings) and duplicate site/model rows.
Fixed-effect inverse-variance pooling assumes a common effect and independent
site estimates. Q, I2, and DerSimonian–Laird random-effects summaries are
descriptive; one-site rows are explicitly labeled and have unavailable
heterogeneity. Site patient counts are summed counts, not deduplicated people
across hospitals. Shared patients across sites would violate independence.
Patient-cluster SEs do not fully address shared calendar/ZCTA exposures; that
and sparse-model calibration require simulation validation.

BH adjustment is performed centrally, with separate families for overall MWAS,
poverty interactions, other SES interactions jointly, and diagnosis-defined
outcomes jointly. The attempted hypothesis count is retained even if some
models cannot be estimated. Reference-SES main effects are not additional
screening hypotheses. `pooled_power_scenarios.csv` sums available null design
information from successful independent-site models; normal power scenarios
are approximations, not BH discovery power or guaranteed sparse-model power.

For a local UCMC demonstration (not multisite evidence):

```bash
Rscript code/28_pool_federated_mwas.R output/mwas/UCMC_first_pass/federated/pool_ucmc output/mwas/UCMC_first_pass/federated/site_estimates.csv
Rscript code/30_audit_duration_candidates.R
Rscript code/29_report_federated_mwas.R
Rscript tests/test_mwas_federated.R
Rscript tests/test_mwas_windows.R
python3 tests/test_acs_ses.py
```

The report accepts `MWAS_RUN_DIR` and `MWAS_POOL_DIR`. The site run can still
produce overall/diagnosis outputs without ACS, with SES explicitly unavailable;
such a protocol cannot be pooled with ACS-enabled outputs.

`code/30_audit_duration_candidates.R` optionally checks single-site overall
FDR-passing candidates against model-based uncertainty and refits after omitting
each patient. It writes aggregate diagnostics only, does not change eligibility
or FDR, and can be expensive for frequent organisms. It cannot validate a
multisite estimate without each site's local diagnostics. The development UCMC
28-day ozone / Candida parapsilosis candidate has only 12 admissions, an extreme
odds ratio, much larger model-based than sandwich uncertainty, and failures in
4/12 leave-one-patient-out fits. It is flagged for review, not a validated signal.
Exploratory clinical-group results use Charlson 0–2 versus 3+, and complete
first-24h modified SOFA 0–5 versus 6+. These cutoffs are development choices;
formal interaction tests are not implemented. SOFA groups define the severity of
organism-positive presentations and are not baseline susceptibility strata.

After models complete, generate the local HTML report:

```bash
MWAS_RUN_DIR=output/mwas/UCMC_first_pass Rscript --vanilla code/24_summarize_acute_mwas.R
```

The report links aggregate tables and plots. `summary.json` provides a machine-
readable overview; the `source/` folder preserves the source used for the run.

## Buddy testing and calibration

See [the buddy guide](../docs/buddy_testing.md) for the frozen protocol and complete installation/run instructions. Scripts 31–36 add count-model uncertainty checks, actual-series simulation calibration, site schema preflight, a resumable local runner, synthetic smoke tests and a source-only bundle. Clinical sensitivities and admission/culture-selection companions are included in script 27. Coordinator script 28 keeps inference methods in separate pooling/FDR families.

Script 37 creates Table 1 and annual site characteristics before script 27. The buddy pilot now requires six core tables, including patient demographics and hospitalization discharge category. Protocol v4 adds seven-day age/sex/race/Charlson/SOFA interactions with preserved coefficient covariance and separate demographic/clinical FDR families. See the buddy guide for references, scales, derivation scopes and missingness.
