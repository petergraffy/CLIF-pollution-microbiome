# Air Pollution, Respiratory Microbial Ecology, and ARF in CLIF

## Working hypothesis

Ambient air pollution alters the culture-detected respiratory microbial ecology of critically ill patients by selecting for some organisms, reducing detection of others, and shaping antimicrobial resistance patterns. These geographically patterned microbial states may increase vulnerability to acute respiratory failure (ARF) and/or ARF severity among patients with pneumonia or sepsis.

Because CLIF contains clinical culture results rather than sequencing-based microbiome assays, this project should be framed as respiratory microbial ecology or culture-detected pathogen community structure, not as a full lung microbiome study.

## Primary aims

1. Estimate whether county-year PM2.5 and NO2 exposure are associated with respiratory culture organism composition among adult ICU hospitalizations.
2. Test whether pollution-associated organism profiles are associated with incident ARF, ARF subtype, and ARF severity.
3. Among pneumonia or sepsis hospitalizations, test whether geographic exposure patterns identify distinct microbial-respiratory failure subphenotypes.

## Core CLIF tables

Required:

- `patient`: demographics and death timestamp.
- `hospitalization`: admission/discharge dates, age, discharge disposition, county/tract/ZIP geography.
- `adt`: first ICU time, ICU length of stay, hospital/site structure.
- `microbiology_culture`: organism, organism group, culture method, specimen fluid, collection time.
- `microbiology_susceptibility`: organism-level antimicrobial susceptibility and resistance phenotypes.
- `respiratory_support`: respiratory support intensity, IMV/NIV/HFNC exposure, ventilator-free day ingredients.

Strongly recommended:

- `labs`: PaO2, PaCO2, pH for physiologic ARF and ARF subtype.
- `vitals`: SpO2 for hypoxemic ARF proxy.
- `hospital_diagnosis`: pneumonia, sepsis, comorbidity, and present-on-admission stratification.
- `medication_admin_continuous`: vasopressor/sepsis severity covariates.

## Exposure linkage

Use `hospitalization.county_code` as the preferred residence county linkage field, with admission year as the exposure year. Existing project exposome files can provide:

- `pm25_county_year.csv` or `conus_county_pm25_2005_2024.csv`
- `no2_county_year.csv` or `conus_county_no2_2005_2024.csv`
- weather covariates from Daymet county-year summaries
- SVI or related county-year covariates where available

Analyses should separately evaluate annual exposure, lagged annual exposure, and, if available, month-of-admission exposure.

## Cohort definition

Base denominator:

- Adult hospitalizations, age >= 18.
- ICU admission identified in `adt`.
- ICU length of stay >= 24 hours for primary analyses.
- Valid CONUS county FIPS.
- Admission year in the available exposure window.

Microbiology analytic cohort:

- At least one respiratory culture around ICU admission or ARF onset.
- Primary culture window: from 48 hours before first ICU admission through 72 hours after first ICU admission.
- Respiratory specimens include `respiratory_tract`, `respiratory_tract_lower`, `nasopharynx_upperairway`, `oropharynx_tongue_oralcavity`, and possibly `pleural_cavity_fluid` for sensitivity analyses.
- Method should preferentially be `culture`; Gram stain and smear can be profiled separately.

Pneumonia/sepsis subcohorts:

- Pneumonia: ICD diagnosis pattern for bacterial/viral/aspiration pneumonia, ideally present on admission when `poa_present` is available.
- Sepsis: ICD sepsis/septic shock codes, plus optional vasopressor/lactate/sepsis-event definitions if available.

## Candidate phenotypes

Microbial features:

- Any positive respiratory culture.
- No growth versus organism detected.
- Organism group indicators.
- Top organism groups by site after minimum-cell suppression.
- Polymicrobial culture indicator.
- Culture-detected diversity metrics such as observed organism groups and Shannon index.
- Antimicrobial non-susceptibility indicators linked by `organism_id`.

ARF and severity outcomes:

- Physiologic ARF using SpO2/FiO2, PaO2/FiO2, PaCO2/pH logic from prior CLIF ARF pollution code.
- Hypoxemic, hypercapnic, and mixed ARF subtypes.
- Highest respiratory support category: none, conventional oxygen, HFNC/NIV, IMV.
- IMV initiation, ICU length of stay, hospital mortality, and ventilator-free days if support intervals are complete enough.

## Initial modeling strategy

1. Site-level QC: organism distributions, specimen fluid mix, culture timing, missing county, exposure coverage.
2. Ecologic microbial models: organism group detection as outcome; PM2.5/NO2 as exposures; adjust for age, sex, race/ethnicity, admission year, season, hospital/site, and specimen type.
3. ARF vulnerability models: ARF or severe ARF as outcome; exposures plus microbial profile features as predictors.
4. Subphenotype models: cluster respiratory culture profiles among pneumonia/sepsis ARF cases, then test exposure gradients and outcome differences.
5. Sensitivity analyses: lower respiratory samples only, early cultures only, positive cultures only, site fixed effects, county random effects, lagged exposures.

## Key limitations to state upfront

- Clinical cultures are ordered based on clinician suspicion, so culture availability and organism detection are confounded by care patterns.
- Culture data under-detect anaerobes, viruses, unculturable organisms, and colonizing flora.
- County-year exposure linkage may misclassify individual exposure and may reflect residence rather than hospital location depending on each site's ETL.
- Antibiotics before culture collection can distort organism detection and susceptibility patterns.

