## Code directory

This directory contains the executable R workflow for the CLIF pollution-microbiome project. The current implementation is exploratory and aggregate-based: it builds a local adult ICU respiratory culture cohort, links county-year PM2.5 and NO2, computes crude organism-pollution correlations, and creates a tile plot.

## Scripts

1. `01_microbiome_cohort_export.R`
   Builds the adult ICU cohort, identifies early respiratory cultures, summarizes organism groups and organism categories, links optional susceptibility data, links county-year exposures, and writes aggregate CSVs.

   Main outputs:
   - `microbe_site_county_year_<site>_<stamp>.csv`
   - `microbe_organism_group_<site>_<stamp>.csv`
   - `microbe_organism_category_<site>_<stamp>.csv`
   - `microbe_qc_summary_<site>_<stamp>.csv`

2. `02_pollution_microbe_correlation.R`
   Reads the latest export files, completes the organism-by-county-year grid so absent organisms are treated as zero counts, and computes Pearson and Spearman correlations for PM2.5 and NO2.

   Main output:
   - `pollution_microbe_correlations_<stamp>.csv`

3. `03_plot_pollution_microbe_correlations.R`
   Reads the latest correlation output and creates a tile plot of organism-category correlations with pollution.

   Main outputs:
   - `pollution_microbe_correlation_tile_<stamp>.png`
   - `pollution_microbe_correlation_tile_<stamp>.pdf`

## Run Order

From the repository root:

```bash
Rscript code/01_microbiome_cohort_export.R
Rscript code/02_pollution_microbe_correlation.R
Rscript code/03_plot_pollution_microbe_correlations.R
```

## Current Analytic Definitions

1. Base cohort: adult ICU hospitalizations with ICU length of stay at least 24 hours and a valid CONUS county FIPS.
2. Culture window: 48 hours before through 72 hours after first ICU admission.
3. Primary respiratory specimens: `respiratory_tract` and `respiratory_tract_lower`.
4. Organism group analysis uses `organism_group`.
5. Genus/species analysis uses `organism_category`.
6. Severity proxy uses early HFNC, NIV/NIPPV/CPAP, or IMV in `respiratory_support`.

## Development Notes

The current code intentionally keeps the exploratory outputs unsuppressed for local signal-finding. Before sharing outside an approved environment, add minimum-cell suppression or review aggregate release rules. The next major code step is to replace the respiratory support proxy with the fuller physiologic ARF phenotype using SpO2/FiO2, PaO2/FiO2, and PaCO2/pH logic.

