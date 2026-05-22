# ================================================================================================
# CLIF Pollution-Microbiome | County-Year Pollution vs Organism Correlations
# Purpose:
#   Use aggregate outputs from 01_microbiome_cohort_export.R to estimate crude county-year
#   correlations between air pollution and organism_category prevalence.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(glue)
})

latest_file <- function(pattern) {
  hits <- Sys.glob(pattern)
  if (length(hits) == 0) stop("No files found for pattern: ", pattern)
  hits[which.max(file.info(hits)$mtime)]
}

safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  unname(cor(x[keep], y[keep], method = method))
}

site_county_year_path <- latest_file("output/final/microbe_site_county_year_*.csv")
organism_category_path <- latest_file("output/final/microbe_organism_category_*.csv")

message("Using site-count county-year file: ", site_county_year_path)
message("Using organism-category file: ", organism_category_path)

site_county_year <- readr::read_csv(site_county_year_path, show_col_types = FALSE) %>%
  mutate(
    county_fips = str_pad(as.character(county_fips), 5, side = "left", pad = "0"),
    year = as.integer(year)
  )

organism_category <- readr::read_csv(organism_category_path, show_col_types = FALSE) %>%
  mutate(
    county_fips = str_pad(as.character(county_fips), 5, side = "left", pad = "0"),
    year = as.integer(year),
    organism_category = as.character(organism_category)
  )

analysis_dat <- tidyr::crossing(
  site_county_year %>%
    select(site_name, county_fips, year, n_icu, n_with_resp_culture,
           n_positive_resp_culture, pm25_mean, no2_mean),
  organism_category = sort(unique(organism_category$organism_category))
) %>%
  left_join(
    organism_category %>%
      select(site_name, county_fips, year, organism_category, n_hospitalizations, n_organism_isolates),
    by = c("site_name", "county_fips", "year", "organism_category")
  ) %>%
  mutate(
    n_hospitalizations = coalesce(n_hospitalizations, 0),
    n_organism_isolates = coalesce(n_organism_isolates, 0),
    prop_among_icu = n_hospitalizations / n_icu,
    prop_among_resp_cultured = n_hospitalizations / n_with_resp_culture,
    prop_among_positive_resp_culture = n_hospitalizations / n_positive_resp_culture
  )

correlation_results <- analysis_dat %>%
  group_by(organism_category) %>%
  summarise(
    total_hospitalizations = sum(n_hospitalizations, na.rm = TRUE),
    n_county_years = n_distinct(paste(county_fips, year)),
    n_county_years_resp_culture_denom = sum(!is.na(prop_among_resp_cultured) & is.finite(prop_among_resp_cultured)),
    n_county_years_positive_denom = sum(!is.na(prop_among_positive_resp_culture) & is.finite(prop_among_positive_resp_culture)),
    pearson_pm25_resp_cultured = safe_cor(pm25_mean, prop_among_resp_cultured, "pearson"),
    spearman_pm25_resp_cultured = safe_cor(pm25_mean, prop_among_resp_cultured, "spearman"),
    pearson_no2_resp_cultured = safe_cor(no2_mean, prop_among_resp_cultured, "pearson"),
    spearman_no2_resp_cultured = safe_cor(no2_mean, prop_among_resp_cultured, "spearman"),
    pearson_pm25_positive = safe_cor(pm25_mean, prop_among_positive_resp_culture, "pearson"),
    spearman_pm25_positive = safe_cor(pm25_mean, prop_among_positive_resp_culture, "spearman"),
    pearson_no2_positive = safe_cor(no2_mean, prop_among_positive_resp_culture, "pearson"),
    spearman_no2_positive = safe_cor(no2_mean, prop_among_positive_resp_culture, "spearman"),
    .groups = "drop"
  ) %>%
  mutate(
    max_abs_spearman_resp_cultured = pmax(
      abs(spearman_pm25_resp_cultured),
      abs(spearman_no2_resp_cultured),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(total_hospitalizations), organism_category)

out_dir <- file.path("output", "final")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- file.path(out_dir, glue("pollution_microbe_correlations_{stamp}.csv"))
readr::write_csv(correlation_results, out_path)

stable_results <- correlation_results %>%
  filter(total_hospitalizations >= 25, n_county_years_resp_culture_denom >= 10)

message("Wrote correlation results: ", out_path)
message("")
message("Top positive Spearman correlations with PM2.5 among respiratory-cultured hospitalizations:")
stable_results %>%
  arrange(desc(spearman_pm25_resp_cultured)) %>%
  select(organism_category, total_hospitalizations, n_county_years_resp_culture_denom,
         spearman_pm25_resp_cultured, spearman_no2_resp_cultured) %>%
  head(12) %>%
  print(n = 12)

message("")
message("Top positive Spearman correlations with NO2 among respiratory-cultured hospitalizations:")
stable_results %>%
  arrange(desc(spearman_no2_resp_cultured)) %>%
  select(organism_category, total_hospitalizations, n_county_years_resp_culture_denom,
         spearman_pm25_resp_cultured, spearman_no2_resp_cultured) %>%
  head(12) %>%
  print(n = 12)

message("")
message("Top negative Spearman correlations with NO2 among respiratory-cultured hospitalizations:")
stable_results %>%
  arrange(spearman_no2_resp_cultured) %>%
  select(organism_category, total_hospitalizations, n_county_years_resp_culture_denom,
         spearman_pm25_resp_cultured, spearman_no2_resp_cultured) %>%
  head(12) %>%
  print(n = 12)
