# ================================================================================================
# CLIF Pollution-Microbiome | Federated Site Export Script
# Purpose:
#   Build an adult ICU cohort, link county-year air pollution exposure files, summarize respiratory
#   culture ecology, and export PHI-safe aggregate tables for pooled analysis.
#
# Required CLIF v2.1+ tables:
#   clif_patient, clif_hospitalization, clif_adt, clif_microbiology_culture,
#   clif_microbiology_susceptibility, clif_respiratory_support
#
# Optional CLIF tables:
#   clif_hospital_diagnosis
#
# Outputs:
#   output/final/microbe_site_county_year_<site>_<stamp>.csv
#   output/final/microbe_organism_group_<site>_<stamp>.csv
#   output/final/microbe_qc_summary_<site>_<stamp>.csv
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(readr)
  library(arrow)
  library(glue)
  library(janitor)
})

source("utils/config.R")

`%||%` <- function(x, y) if (!is.null(x)) x else y

site_name <- config$site_name
tables_path <- config$tables_path
file_type <- tolower(config$file_type)
exposome_path <- config$exposome_path %||% file.path(config$repo %||% getwd(), "exposome")

START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
END_DATE <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")
MIN_ICU_LOS_H <- 24
MICROBE_WINDOW_PRE_H <- 48
MICROBE_WINDOW_POST_H <- 72
RESP_FLUID_PRIMARY <- c("respiratory_tract", "respiratory_tract_lower")
RESP_FLUID_SENSITIVITY <- c(
  RESP_FLUID_PRIMARY,
  "nasopharynx_upperairway",
  "oropharynx_tongue_oralcavity",
  "pleural_cavity_fluid"
)

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

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

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(
    ext,
    "csv" = readr::read_csv(path, show_col_types = FALSE),
    "parquet" = arrow::read_parquet(path),
    "fst" = {
      if (!requireNamespace("fst", quietly = TRUE)) {
        stop("Package 'fst' is required to read fst files.")
      }
      fst::read_fst(path)
    },
    stop("Unsupported extension: ", ext, " for path: ", path)
  ) %>%
    janitor::clean_names()
}

find_table_path <- function(tbl_base, required = TRUE) {
  wanted <- tolower(tbl_base)
  if (!startsWith(wanted, "clif_")) wanted <- paste0("clif_", wanted)
  files <- list.files(tables_path, full.names = TRUE, recursive = TRUE)
  files <- files[grepl("\\.(csv|parquet|fst)$", files, ignore.case = TRUE)]
  base <- tolower(tools::file_path_sans_ext(basename(files)))
  base_norm <- ifelse(startsWith(base, "clif_"), base, paste0("clif_", base))
  hit <- files[base_norm == wanted]
  if (length(hit) == 1) return(hit)
  if (required) stop(glue("Could not uniquely locate {wanted} in {tables_path}. Matches: {length(hit)}"))
  NA_character_
}

get_tbl <- function(tbl_base, required = TRUE) {
  path <- find_table_path(tbl_base, required = required)
  if (is.na(path)) return(NULL)
  read_any(path)
}

normalize_county_fips <- function(x) {
  x <- str_replace_all(as.character(x), "[^0-9]", "")
  x <- ifelse(nchar(x) > 0 & nchar(x) <= 5, str_pad(x, 5, side = "left", pad = "0"), x)
  ifelse(nchar(x) == 5, x, NA_character_)
}

is_conus_county_fips <- function(x) {
  x <- normalize_county_fips(x)
  state_fips <- substr(x, 1, 2)
  !is.na(x) & !(state_fips %in% c("02", "15", "60", "66", "69", "72", "78"))
}

age_band_4 <- function(age) {
  cut(
    age,
    breaks = c(18, 40, 65, 75, Inf),
    right = FALSE,
    labels = c("18-39", "40-64", "65-74", "75+")
  )
}

harmonize_sex <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x %in% c("female", "f") ~ "Female",
    x %in% c("male", "m") ~ "Male",
    TRUE ~ "Other/Unknown"
  )
}

load_pollutant <- function(pollutant) {
  candidates <- switch(
    pollutant,
    "pm25" = c("conus_county_pm25_2005_2024.csv", "pm25_county_year.csv"),
    "no2" = c("conus_county_no2_2005_2024.csv", "no2_county_year.csv")
  )
  path <- file.path(exposome_path, candidates)
  path <- path[file.exists(path)][1]
  if (is.na(path)) {
    warning("No ", pollutant, " county-year exposure file found in ", exposome_path)
    out <- tibble(county_fips = character(), year = integer(), value = numeric())
    return(rename(out, !!paste0(pollutant, "_mean") := value))
  }
  df <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  value_col <- switch(
    pollutant,
    "pm25" = names(df)[names(df) %in% c("pm25_mean", "mean_pm25", "pm2_5_mean", "pm25", "mean_pm2_5")][1],
    "no2" = names(df)[names(df) %in% c("no2_mean", "mean_no2", "no2", "mean_no_2", "avg_no2")][1]
  )
  geoid_col <- names(df)[names(df) %in% c("geoid", "county_fips", "fips", "county_code")][1]
  year_col <- names(df)[names(df) %in% c("year", "yr")][1]
  if (is.na(value_col) || is.na(geoid_col) || is.na(year_col)) {
    stop("Exposure file is missing geoid/year/value columns: ", path)
  }
  df %>%
    transmute(
      county_fips = normalize_county_fips(.data[[geoid_col]]),
      year = as.integer(.data[[year_col]]),
      value = as.numeric(.data[[value_col]])
    ) %>%
    distinct(county_fips, year, .keep_all = TRUE) %>%
    rename(!!paste0(pollutant, "_mean") := value)
}

message(glue("Site: {site_name}"))
message(glue("Tables path: {tables_path}"))
message(glue("Exposure path: {exposome_path}"))

patient <- get_tbl("patient") %>%
  transmute(
    patient_id,
    sex_category,
    race_category,
    ethnicity_category,
    death_dttm = safe_ts(death_dttm)
  )

hospitalization <- get_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission)),
    discharge_category,
    county_fips = normalize_county_fips(county_code)
  )

adt <- get_tbl("adt") %>%
  transmute(
    hospitalization_id,
    hospital_id = as.character(hospital_id),
    in_dttm = safe_ts(in_dttm),
    out_dttm = safe_ts(out_dttm),
    location_category = str_to_lower(as.character(location_category))
  )

icu_bounds <- adt %>%
  filter(location_category == "icu", !is.na(in_dttm), !is.na(out_dttm), out_dttm > in_dttm) %>%
  group_by(hospitalization_id) %>%
  summarise(
    first_icu_in = min(in_dttm),
    last_icu_out = max(out_dttm),
    hospital_id = first(hospital_id[!is.na(hospital_id)], default = NA_character_),
    icu_los_hours = as.numeric(difftime(last_icu_out, first_icu_in, units = "hours")),
    .groups = "drop"
  )

base <- hospitalization %>%
  inner_join(icu_bounds, by = "hospitalization_id") %>%
  left_join(patient, by = "patient_id") %>%
  mutate(
    year = year(admission_dttm),
    age_band = age_band_4(age_at_admission),
    sex = harmonize_sex(sex_category),
    adult = !is.na(age_at_admission) & age_at_admission >= 18,
    icu_24h = !is.na(icu_los_hours) & icu_los_hours >= MIN_ICU_LOS_H,
    conus_county = is_conus_county_fips(county_fips),
    include_base = adult & icu_24h & conus_county &
      first_icu_in >= START_DATE & first_icu_in <= END_DATE
  )

base_ids <- base %>% filter(include_base) %>% pull(hospitalization_id)

micro <- get_tbl("microbiology_culture") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    patient_id,
    hospitalization_id,
    organism_id = as.character(organism_id),
    collect_dttm = safe_ts(collect_dttm),
    fluid_category = str_to_lower(str_trim(as.character(fluid_category))),
    method_category = str_to_lower(str_trim(as.character(method_category))),
    organism_category = str_to_lower(str_trim(as.character(organism_category))),
    organism_group = str_to_lower(str_trim(as.character(organism_group)))
  ) %>%
  mutate(
    organism_group = na_if(organism_group, ""),
    organism_group = coalesce(organism_group, organism_category),
    no_growth = organism_group %in% c("no_growth", "no growth"),
    positive_culture = !is.na(organism_group) & !no_growth
  )

micro_win <- base %>%
  filter(include_base) %>%
  transmute(
    hospitalization_id,
    culture_window_start = first_icu_in - hours(MICROBE_WINDOW_PRE_H),
    culture_window_end = first_icu_in + hours(MICROBE_WINDOW_POST_H)
  ) %>%
  inner_join(micro, by = "hospitalization_id") %>%
  filter(!is.na(collect_dttm), collect_dttm >= culture_window_start, collect_dttm <= culture_window_end) %>%
  mutate(
    primary_resp_sample = fluid_category %in% RESP_FLUID_PRIMARY,
    sensitivity_resp_sample = fluid_category %in% RESP_FLUID_SENSITIVITY,
    culture_method = method_category == "culture"
  )

susceptibility <- get_tbl("microbiology_susceptibility", required = FALSE)
if (!is.null(susceptibility)) {
  susceptibility <- susceptibility %>%
    transmute(
      organism_id = as.character(organism_id),
      antimicrobial_category = str_to_lower(str_trim(as.character(antimicrobial_category))),
      susceptibility_category = str_to_lower(str_trim(as.character(susceptibility_category)))
    )

  organism_resistance <- susceptibility %>%
    group_by(organism_id) %>%
    summarise(
      n_antimicrobials_tested = n_distinct(antimicrobial_category, na.rm = TRUE),
      any_non_susceptible = any(susceptibility_category == "non_susceptible", na.rm = TRUE),
      .groups = "drop"
    )
} else {
  warning("microbiology_susceptibility table not found; resistance exports will be unavailable.")
  organism_resistance <- tibble(
    organism_id = character(),
    n_antimicrobials_tested = integer(),
    any_non_susceptible = logical()
  )
}

micro_features <- micro_win %>%
  filter(primary_resp_sample, culture_method) %>%
  left_join(organism_resistance, by = "organism_id") %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_resp_culture = TRUE,
    any_positive_resp_culture = any(positive_culture, na.rm = TRUE),
    any_no_growth = any(no_growth, na.rm = TRUE),
    n_resp_cultures = n_distinct(collect_dttm, fluid_category, method_category, na.rm = TRUE),
    n_positive_organism_groups = n_distinct(organism_group[positive_culture], na.rm = TRUE),
    polymicrobial = n_positive_organism_groups >= 2,
    any_non_susceptible = any(any_non_susceptible, na.rm = TRUE),
    .groups = "drop"
  )

resp_support <- get_tbl("respiratory_support") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    device_category = str_to_lower(str_trim(as.character(device_category)))
  )

resp_features <- resp_support %>%
  inner_join(
    base %>% filter(include_base) %>% select(hospitalization_id, first_icu_in),
    by = "hospitalization_id"
  ) %>%
  filter(recorded_dttm >= first_icu_in - hours(24), recorded_dttm <= first_icu_in + hours(72)) %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_hfnc = any(device_category %in% c("hfnc", "high flow nc", "high_flow_nasal_cannula"), na.rm = TRUE),
    any_niv = any(device_category %in% c("niv", "nippv", "bipap", "cpap"), na.rm = TRUE),
    any_imv = any(device_category %in% c("imv", "invasive_mechanical_ventilation"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    severe_arf_proxy = any_imv | any_niv | any_hfnc,
    max_support = case_when(
      any_imv ~ "imv",
      any_niv ~ "niv",
      any_hfnc ~ "hfnc",
      TRUE ~ "none_or_conventional"
    )
  )

diagnosis <- get_tbl("hospital_diagnosis", required = FALSE)
if (!is.null(diagnosis)) {
  diagnosis_flags <- diagnosis %>%
    filter(hospitalization_id %in% base_ids) %>%
    mutate(
      diagnosis_code = str_to_upper(str_replace_all(as.character(diagnosis_code), "\\.", "")),
      pneumonia_dx = str_detect(diagnosis_code, "^(J12|J13|J14|J15|J16|J17|J18|J69)"),
      sepsis_dx = str_detect(diagnosis_code, "^(A40|A41|R652|R65)")
    ) %>%
    group_by(hospitalization_id) %>%
    summarise(
      pneumonia_dx = any(pneumonia_dx, na.rm = TRUE),
      sepsis_dx = any(sepsis_dx, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  diagnosis_flags <- tibble(
    hospitalization_id = base_ids,
    pneumonia_dx = NA,
    sepsis_dx = NA
  )
}

exposure <- load_pollutant("pm25") %>%
  full_join(load_pollutant("no2"), by = c("county_fips", "year"))

analysis_cohort <- base %>%
  filter(include_base) %>%
  left_join(micro_features, by = "hospitalization_id") %>%
  left_join(resp_features, by = "hospitalization_id") %>%
  left_join(diagnosis_flags, by = "hospitalization_id") %>%
  left_join(exposure, by = c("county_fips", "year")) %>%
  mutate(
    across(
      c(any_resp_culture, any_positive_resp_culture, any_no_growth, polymicrobial,
        any_non_susceptible, any_hfnc, any_niv, any_imv, severe_arf_proxy),
      ~ coalesce(.x, FALSE)
    ),
    max_support = coalesce(max_support, "none_or_conventional"),
    in_hospital_death = discharge_category == "Expired" | (!is.na(death_dttm) & death_dttm <= discharge_dttm)
  )

site_county_year <- analysis_cohort %>%
  group_by(site_name = site_name, county_fips, year) %>%
  summarise(
    n_icu = n(),
    n_with_resp_culture = sum(any_resp_culture),
    n_positive_resp_culture = sum(any_positive_resp_culture),
    n_polymicrobial = sum(polymicrobial),
    n_non_susceptible = if (!is.null(susceptibility)) sum(any_non_susceptible) else NA_integer_,
    n_severe_arf_proxy = sum(severe_arf_proxy),
    n_imv = sum(any_imv),
    n_deaths = sum(in_hospital_death, na.rm = TRUE),
    pm25_mean = first(pm25_mean),
    no2_mean = first(no2_mean),
    .groups = "drop"
  )

organism_group_summary <- micro_win %>%
  filter(primary_resp_sample, culture_method, positive_culture) %>%
  inner_join(analysis_cohort %>% select(hospitalization_id, county_fips, year), by = "hospitalization_id") %>%
  group_by(site_name = site_name, county_fips, year, organism_group) %>%
  summarise(
    n_hospitalizations = n_distinct(hospitalization_id),
    n_organism_isolates = n(),
    .groups = "drop"
  )

organism_category_summary <- micro_win %>%
  filter(primary_resp_sample, culture_method, positive_culture) %>%
  inner_join(analysis_cohort %>% select(hospitalization_id, county_fips, year), by = "hospitalization_id") %>%
  group_by(site_name = site_name, county_fips, year, organism_category) %>%
  summarise(
    n_hospitalizations = n_distinct(hospitalization_id),
    n_organism_isolates = n(),
    .groups = "drop"
  )

qc_summary <- tibble(
  site_name = site_name,
  n_hospitalizations = nrow(hospitalization),
  n_adult_icu_24h_conus = nrow(analysis_cohort),
  n_missing_county = sum(is.na(base$county_fips)),
  n_with_primary_resp_culture = sum(analysis_cohort$any_resp_culture),
  n_with_positive_primary_resp_culture = sum(analysis_cohort$any_positive_resp_culture),
  susceptibility_available = !is.null(susceptibility),
  n_with_pneumonia_dx = sum(analysis_cohort$pneumonia_dx, na.rm = TRUE),
  n_with_sepsis_dx = sum(analysis_cohort$sepsis_dx, na.rm = TRUE),
  n_with_pm25 = sum(!is.na(analysis_cohort$pm25_mean)),
  n_with_no2 = sum(!is.na(analysis_cohort$no2_mean))
)

readr::write_csv(
  site_county_year,
  file.path(out_dir, glue("microbe_site_county_year_{site_name}_{stamp}.csv"))
)
readr::write_csv(
  organism_group_summary,
  file.path(out_dir, glue("microbe_organism_group_{site_name}_{stamp}.csv"))
)
readr::write_csv(
  organism_category_summary,
  file.path(out_dir, glue("microbe_organism_category_{site_name}_{stamp}.csv"))
)
readr::write_csv(
  qc_summary,
  file.path(out_dir, glue("microbe_qc_summary_{site_name}_{stamp}.csv"))
)

message(glue("Wrote PHI-safe aggregate exports to {out_dir}"))
