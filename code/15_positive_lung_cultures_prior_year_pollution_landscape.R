# ================================================================================================
# Positive Lung Culture Prior-Year ZCTA Pollution Landscape
#
# Goal:
#   Export every positive pulmonary culture row in the CLIF database and link each hospitalization
#   to ZCTA-level PM2.5 and NO2 in the calendar year before hospital admission.
#
# Notes:
#   - This is a landscape/export script, not a modeling script.
#   - "Lung" uses the project's primary pulmonary fluid categories:
#       respiratory_tract and respiratory_tract_lower.
#   - Positive culture means method_category == culture and organism_group/category is not no_growth.
#   - Prior-year exposure is admission_year - 1.
# ================================================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(glue)
  library(janitor)
  library(lubridate)
  library(readr)
  library(stringr)
})

source("utils/clif_io.R")

site_name <- clif_site_name
tables_path <- clif_tables_path
zcta_dir <- clif_zcta_exposure_dir

pm25_path <- find_zcta_exposure_path("air_pollution_zcta_pm25_monthly_2005_2023.parquet")
no2_path <- find_zcta_exposure_path("air_pollution_zcta_no2_annual_2005_2025.parquet")

safe_ts <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = tz))
  if (is.numeric(x)) {
    x2 <- ifelse(x > 1e12, x / 1000, x)
    return(as.POSIXct(x2, origin = "1970-01-01", tz = tz))
  }
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c("ymd_HMS", "ymd_HM", "ymd", "ymdTz", "ymdT", "mdy_HMS", "mdy_HM", "mdy"),
    tz = tz,
    quiet = TRUE
  ))
}

normalize_zip <- function(x) {
  x <- str_replace_all(as.character(x), "[^0-9]", "")
  x <- ifelse(nchar(x) >= 5, substr(x, 1, 5), x)
  ifelse(nchar(x) == 5, x, NA_character_)
}

count_nonmissing <- function(x) sum(!is.na(x))

message("Using CLIF tables: ", tables_path)
message("Using ZCTA exposures: ", zcta_dir)

pulmonary_primary <- c("respiratory_tract", "respiratory_tract_lower")

patient <- read_tbl("patient") %>%
  transmute(
    patient_id,
    sex_category = str_to_lower(str_trim(as.character(sex_category))),
    race_category = str_to_lower(str_trim(as.character(race_category))),
    ethnicity_category = str_to_lower(str_trim(as.character(ethnicity_category)))
  )

hospitalization <- read_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm),
    admission_year = year(admission_dttm),
    exposure_year = admission_year - 1L,
    age_at_admission = suppressWarnings(as.numeric(age_at_admission)),
    zipcode_five_digit = normalize_zip(zipcode_five_digit),
    state_code = str_to_upper(str_trim(as.character(state_code))),
    county_code = str_trim(as.character(county_code))
  ) %>%
  left_join(patient, by = "patient_id")

positive_lung_cultures <- read_tbl("microbiology_culture") %>%
  transmute(
    patient_id,
    hospitalization_id,
    organism_id,
    order_dttm = safe_ts(order_dttm),
    collect_dttm = safe_ts(collect_dttm),
    result_dttm = safe_ts(result_dttm),
    fluid_name = str_to_lower(str_trim(as.character(fluid_name))),
    fluid_category = str_to_lower(str_trim(as.character(fluid_category))),
    method_name = str_to_lower(str_trim(as.character(method_name))),
    method_category = str_to_lower(str_trim(as.character(method_category))),
    organism_name = str_to_lower(str_trim(as.character(organism_name))),
    organism_category = str_to_lower(str_trim(as.character(organism_category))),
    organism_group = str_to_lower(str_trim(as.character(organism_group)))
  ) %>%
  mutate(
    organism_group = coalesce(na_if(organism_group, ""), organism_category),
    no_growth = organism_group %in% c("no_growth", "no growth"),
    positive_culture = !is.na(organism_group) & !no_growth,
    pulmonary_primary = fluid_category %in% pulmonary_primary
  ) %>%
  filter(method_category == "culture", pulmonary_primary, positive_culture) %>%
  mutate(culture_row_id = row_number(), .before = patient_id)

pm25_prior_year <- arrow::read_parquet(pm25_path) %>%
  transmute(
    zipcode_five_digit = normalize_zip(zip),
    exposure_year = as.integer(year),
    pm25_prior_year = as.numeric(pm25_ug_m3)
  ) %>%
  group_by(zipcode_five_digit, exposure_year) %>%
  summarise(pm25_prior_year = mean(pm25_prior_year, na.rm = TRUE), .groups = "drop") %>%
  mutate(pm25_prior_year = if_else(is.nan(pm25_prior_year), NA_real_, pm25_prior_year))

no2_prior_year <- arrow::read_parquet(no2_path) %>%
  transmute(
    zipcode_five_digit = normalize_zip(zip),
    exposure_year = as.integer(year),
    no2_prior_year = as.numeric(no2)
  ) %>%
  group_by(zipcode_five_digit, exposure_year) %>%
  summarise(no2_prior_year = mean(no2_prior_year, na.rm = TRUE), .groups = "drop") %>%
  mutate(no2_prior_year = if_else(is.nan(no2_prior_year), NA_real_, no2_prior_year))

culture_exposure <- positive_lung_cultures %>%
  left_join(hospitalization, by = c("patient_id", "hospitalization_id")) %>%
  left_join(pm25_prior_year, by = c("zipcode_five_digit", "exposure_year")) %>%
  left_join(no2_prior_year, by = c("zipcode_five_digit", "exposure_year")) %>%
  arrange(admission_dttm, hospitalization_id, collect_dttm, organism_category, organism_group)

coverage_summary <- tibble(
  site_name = site_name,
  n_positive_lung_culture_rows = nrow(culture_exposure),
  n_hospitalizations = n_distinct(culture_exposure$hospitalization_id),
  n_patients = n_distinct(culture_exposure$patient_id),
  n_organism_categories = n_distinct(culture_exposure$organism_category, na.rm = TRUE),
  n_organism_groups = n_distinct(culture_exposure$organism_group, na.rm = TRUE),
  n_with_admission_dttm = count_nonmissing(culture_exposure$admission_dttm),
  n_with_zip = count_nonmissing(culture_exposure$zipcode_five_digit),
  n_with_pm25_prior_year = count_nonmissing(culture_exposure$pm25_prior_year),
  n_with_no2_prior_year = count_nonmissing(culture_exposure$no2_prior_year),
  n_with_both_prior_year_exposures = sum(!is.na(culture_exposure$pm25_prior_year) & !is.na(culture_exposure$no2_prior_year)),
  min_admission_year = suppressWarnings(min(culture_exposure$admission_year, na.rm = TRUE)),
  max_admission_year = suppressWarnings(max(culture_exposure$admission_year, na.rm = TRUE)),
  min_exposure_year = suppressWarnings(min(culture_exposure$exposure_year, na.rm = TRUE)),
  max_exposure_year = suppressWarnings(max(culture_exposure$exposure_year, na.rm = TRUE))
) %>%
  mutate(across(starts_with("min_") | starts_with("max_"), ~ if_else(is.infinite(.x), NA_real_, as.numeric(.x))))

organism_summary <- culture_exposure %>%
  distinct(hospitalization_id, organism_category, organism_group, .keep_all = TRUE) %>%
  group_by(organism_category, organism_group) %>%
  summarise(
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    n_culture_rows = n(),
    n_zctas = n_distinct(zipcode_five_digit, na.rm = TRUE),
    n_with_pm25_prior_year = count_nonmissing(pm25_prior_year),
    pm25_prior_year_median = median(pm25_prior_year, na.rm = TRUE),
    pm25_prior_year_iqr = IQR(pm25_prior_year, na.rm = TRUE),
    n_with_no2_prior_year = count_nonmissing(no2_prior_year),
    no2_prior_year_median = median(no2_prior_year, na.rm = TRUE),
    no2_prior_year_iqr = IQR(no2_prior_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pm25_prior_year_median = if_else(is.nan(pm25_prior_year_median), NA_real_, pm25_prior_year_median),
    pm25_prior_year_iqr = if_else(is.nan(pm25_prior_year_iqr), NA_real_, pm25_prior_year_iqr),
    no2_prior_year_median = if_else(is.nan(no2_prior_year_median), NA_real_, no2_prior_year_median),
    no2_prior_year_iqr = if_else(is.nan(no2_prior_year_iqr), NA_real_, no2_prior_year_iqr)
  ) %>%
  arrange(desc(n_hospitalizations), organism_category, organism_group)

year_summary <- culture_exposure %>%
  distinct(hospitalization_id, admission_year, exposure_year, zipcode_five_digit, pm25_prior_year, no2_prior_year) %>%
  group_by(admission_year, exposure_year) %>%
  summarise(
    n_hospitalizations = n_distinct(hospitalization_id),
    n_zctas = n_distinct(zipcode_five_digit, na.rm = TRUE),
    n_with_pm25_prior_year = count_nonmissing(pm25_prior_year),
    pm25_prior_year_median = median(pm25_prior_year, na.rm = TRUE),
    n_with_no2_prior_year = count_nonmissing(no2_prior_year),
    no2_prior_year_median = median(no2_prior_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pm25_prior_year_median = if_else(is.nan(pm25_prior_year_median), NA_real_, pm25_prior_year_median),
    no2_prior_year_median = if_else(is.nan(no2_prior_year_median), NA_real_, no2_prior_year_median)
  ) %>%
  arrange(admission_year)

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

culture_path <- file.path(out_dir, glue("positive_lung_cultures_prior_year_pollution_{site_name}_{stamp}.csv"))
coverage_path <- file.path(out_dir, glue("positive_lung_cultures_prior_year_pollution_coverage_{site_name}_{stamp}.csv"))
organism_path <- file.path(out_dir, glue("positive_lung_cultures_prior_year_pollution_organism_summary_{site_name}_{stamp}.csv"))
year_path <- file.path(out_dir, glue("positive_lung_cultures_prior_year_pollution_year_summary_{site_name}_{stamp}.csv"))

readr::write_csv(culture_exposure, culture_path)
readr::write_csv(coverage_summary, coverage_path)
readr::write_csv(organism_summary, organism_path)
readr::write_csv(year_summary, year_path)

message("Coverage:")
print(coverage_summary)
message("")
message("Top organisms by hospitalization count:")
print(organism_summary %>% select(organism_category, organism_group, n_hospitalizations, n_patients, n_with_pm25_prior_year, n_with_no2_prior_year) %>% head(25), n = 25)
message("")
message("Wrote culture-level export: ", culture_path)
message("Wrote coverage summary: ", coverage_path)
message("Wrote organism summary: ", organism_path)
message("Wrote year summary: ", year_path)
