# ================================================================================================
# SHRF ZCTA Pollution vs Positive Pulmonary Culture Models
#
# Cohort:
#   Adult ED -> ICU hospitalizations with severe hypoxemic respiratory failure:
#   IMV + PaO2/FiO2 < 300 in first 24h after first ICU admission.
#   ED intubations are eligible.
#
# Exposure:
#   ZCTA-level air pollution release:
#   https://github.com/petergraffy/environment_transplant_survival/releases/tag/air-pollution-zcta-v1
#   - PM2.5 monthly 2005-2023, annualized to admission year
#   - Ozone monthly 2005-2023, annualized to admission year
#   - NO2 annual 2005-2025
#
# Outcome:
#   Any positive pulmonary microbiology culture, primary definition:
#   fluid_category in respiratory_tract or respiratory_tract_lower.
#
# Model:
#   outcome ~ pollutant per IQR + age_band + sex + admission_year + (1 | zipcode_five_digit)
# ================================================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(forcats)
  library(glmmTMB)
  library(glue)
  library(janitor)
  library(lubridate)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("utils/clif_io.R")

site_name <- clif_site_name
tables_path <- clif_tables_path
zcta_dir <- clif_zcta_exposure_dir
pm25_path <- find_zcta_exposure_path("air_pollution_zcta_pm25_monthly_2005_2023.parquet")
o3_path <- find_zcta_exposure_path("air_pollution_zcta_o3_monthly_2005_2023.parquet")
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

normalize_fio2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  case_when(
    x > 1 & x <= 100 ~ x / 100,
    x >= 0.21 & x <= 1 ~ x,
    TRUE ~ NA_real_
  )
}

normalize_zip <- function(x) {
  x <- str_replace_all(as.character(x), "[^0-9]", "")
  x <- ifelse(nchar(x) >= 5, substr(x, 1, 5), x)
  ifelse(nchar(x) == 5, x, NA_character_)
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

fit_pollution_model <- function(data, outcome, exposure) {
  dat <- data %>%
    mutate(
      outcome = as.integer(.data[[outcome]]),
      exposure_value = .data[[exposure]]
    ) %>%
    filter(
      !is.na(outcome),
      !is.na(exposure_value),
      !is.na(zipcode_five_digit),
      !is.na(admission_year)
    )

  exposure_iqr <- IQR(dat$exposure_value, na.rm = TRUE)
  if (!is.finite(exposure_iqr) || exposure_iqr <= 0) return(NULL)

  dat <- dat %>%
    mutate(
      exposure_iqr_scaled = exposure_value / exposure_iqr,
      admission_year = factor(admission_year),
      zipcode_five_digit = factor(zipcode_five_digit),
      age_band = fct_na_value_to_level(factor(age_band), level = "Unknown"),
      sex = factor(sex)
    )

  fit <- glmmTMB(
    outcome ~ exposure_iqr_scaled + age_band + sex + admission_year + (1 | zipcode_five_digit),
    family = binomial(),
    data = dat,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )

  coef_tab <- summary(fit)$coefficients$cond
  beta <- coef_tab["exposure_iqr_scaled", "Estimate"]
  se <- coef_tab["exposure_iqr_scaled", "Std. Error"]
  p_value <- coef_tab["exposure_iqr_scaled", "Pr(>|z|)"]

  tibble(
    outcome = outcome,
    exposure = exposure,
    n_hospitalizations = nrow(dat),
    n_events = sum(dat$outcome == 1, na.rm = TRUE),
    n_zctas = n_distinct(dat$zipcode_five_digit),
    exposure_iqr = exposure_iqr,
    odds_ratio_per_iqr = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    p_value = p_value
  )
}

message("Using CLIF tables: ", tables_path)
message("Using ZCTA exposures: ", zcta_dir)

patient <- read_tbl("patient") %>%
  transmute(patient_id, sex_category)

hospitalization <- read_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm),
    admission_year = year(admission_dttm),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission)),
    zipcode_five_digit = normalize_zip(zipcode_five_digit)
  ) %>%
  left_join(patient, by = "patient_id") %>%
  mutate(
    age_band = age_band_4(age_at_admission),
    sex = harmonize_sex(sex_category)
  )

adt <- read_tbl("adt") %>%
  transmute(
    hospitalization_id,
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
    .groups = "drop"
  )

pathway_flags <- adt %>%
  inner_join(icu_bounds %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
  filter(location_category != "icu", !is.na(out_dttm), out_dttm <= first_icu_in) %>%
  group_by(hospitalization_id) %>%
  slice_max(out_dttm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    hospitalization_id,
    pre_icu_location_category = location_category,
    pre_icu_out_dttm = out_dttm,
    ed_to_icu = pre_icu_location_category == "ed" & pre_icu_out_dttm == first_icu_in
  )

base <- hospitalization %>%
  inner_join(icu_bounds, by = "hospitalization_id") %>%
  left_join(pathway_flags, by = "hospitalization_id") %>%
  mutate(
    adult = !is.na(age_at_admission) & age_at_admission >= 18,
    ed_to_icu = coalesce(ed_to_icu, FALSE),
    win_24_start = first_icu_in,
    win_24_end = first_icu_in + hours(24),
    win_48_start = first_icu_in,
    win_48_end = first_icu_in + hours(48),
    early_window_start = first_icu_in - hours(48),
    early_window_end = first_icu_in + hours(72)
  ) %>%
  filter(adult)

base_ids <- base$hospitalization_id

resp_support_all <- read_tbl("respiratory_support") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    device_category = str_to_lower(str_trim(as.character(device_category))),
    fio2_set = normalize_fio2(fio2_set)
  )

pre_icu_imv_flags <- resp_support_all %>%
  inner_join(base %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_before_icu = any(
      device_category %in% c("imv", "invasive_mechanical_ventilation") & recorded_dttm < first_icu_in,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

analysis_base <- base %>%
  left_join(pre_icu_imv_flags, by = "hospitalization_id") %>%
  mutate(imv_before_icu = coalesce(imv_before_icu, FALSE)) %>%
  filter(ed_to_icu)

analysis_base_ids <- analysis_base$hospitalization_id

resp_support <- resp_support_all %>%
  filter(hospitalization_id %in% analysis_base_ids) %>%
  inner_join(analysis_base %>% select(hospitalization_id, win_24_start, win_24_end), by = "hospitalization_id") %>%
  filter(recorded_dttm >= win_24_start, recorded_dttm <= win_24_end)

imv_flags <- resp_support %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_imv_24h = any(device_category %in% c("imv", "invasive_mechanical_ventilation"), na.rm = TRUE),
    .groups = "drop"
  )

labs <- read_tbl("labs") %>%
  filter(hospitalization_id %in% analysis_base_ids) %>%
  transmute(
    hospitalization_id,
    lab_result_dttm = safe_ts(lab_result_dttm),
    lab_category = str_to_lower(as.character(lab_category)),
    lab_value_numeric = suppressWarnings(as.numeric(lab_value_numeric))
  ) %>%
  filter(lab_category == "po2_arterial", !is.na(lab_result_dttm), !is.na(lab_value_numeric)) %>%
  inner_join(analysis_base %>% select(hospitalization_id, win_24_start, win_24_end), by = "hospitalization_id") %>%
  filter(lab_result_dttm >= win_24_start, lab_result_dttm <= win_24_end)

po2_dt <- as.data.table(labs %>% transmute(hospitalization_id, po2_time = lab_result_dttm, po2 = lab_value_numeric))
fio2_dt <- as.data.table(resp_support %>% filter(!is.na(fio2_set)) %>% transmute(hospitalization_id, fio2_time = recorded_dttm, fio2_set))

setkey(po2_dt, hospitalization_id, po2_time)
setkey(fio2_dt, hospitalization_id, fio2_time)
po2_dt[, po2_time_keep := po2_time]

pf_pairs <- fio2_dt[
  po2_dt,
  roll = "nearest",
  on = .(hospitalization_id, fio2_time = po2_time),
  nomatch = 0L
][
  , time_diff_h := abs(as.numeric(difftime(po2_time_keep, fio2_time, units = "hours")))
][
  time_diff_h <= 1,
  .(hospitalization_id, pf_ratio = po2 / fio2_set)
]

pf_flags <- as_tibble(pf_pairs) %>%
  group_by(hospitalization_id) %>%
  summarise(any_pf_lt_300_24h = any(pf_ratio < 300, na.rm = TRUE), .groups = "drop")

shrf_cohort <- analysis_base %>%
  left_join(imv_flags, by = "hospitalization_id") %>%
  left_join(pf_flags, by = "hospitalization_id") %>%
  mutate(
    any_imv_24h = coalesce(any_imv_24h, FALSE),
    any_pf_lt_300_24h = coalesce(any_pf_lt_300_24h, FALSE),
    severe_hypoxemic_rf_24h = any_imv_24h & any_pf_lt_300_24h
  ) %>%
  filter(severe_hypoxemic_rf_24h)

pulmonary_primary <- c("respiratory_tract", "respiratory_tract_lower")

micro <- read_tbl("microbiology_culture") %>%
  filter(hospitalization_id %in% shrf_cohort$hospitalization_id) %>%
  transmute(
    patient_id,
    hospitalization_id,
    collect_dttm = safe_ts(collect_dttm),
    fluid_category = str_to_lower(str_trim(as.character(fluid_category))),
    method_category = str_to_lower(str_trim(as.character(method_category))),
    organism_category = str_to_lower(str_trim(as.character(organism_category))),
    organism_group = str_to_lower(str_trim(as.character(organism_group)))
  ) %>%
  mutate(
    organism_group = coalesce(na_if(organism_group, ""), organism_category),
    no_growth = organism_group %in% c("no_growth", "no growth"),
    positive_culture = !is.na(organism_group) & !no_growth,
    pulmonary_primary = fluid_category %in% pulmonary_primary
  ) %>%
  filter(method_category == "culture", positive_culture, pulmonary_primary)

culture_flags <- shrf_cohort %>%
  select(patient_id, hospitalization_id, admission_dttm, discharge_dttm, first_icu_in,
         win_24_start, win_24_end, win_48_start, win_48_end,
         early_window_start, early_window_end) %>%
  left_join(
    micro %>%
      group_by(hospitalization_id) %>%
      summarise(any_positive_pulm_culture_hosp = TRUE, .groups = "drop"),
    by = "hospitalization_id"
  ) %>%
  left_join(
    micro %>%
      inner_join(shrf_cohort %>% select(hospitalization_id, win_24_start, win_24_end), by = "hospitalization_id") %>%
      filter(collect_dttm >= win_24_start, collect_dttm <= win_24_end) %>%
      group_by(hospitalization_id) %>%
      summarise(any_positive_pulm_culture_24h = TRUE, .groups = "drop"),
    by = "hospitalization_id"
  ) %>%
  left_join(
    micro %>%
      inner_join(shrf_cohort %>% select(hospitalization_id, win_48_start, win_48_end), by = "hospitalization_id") %>%
      filter(collect_dttm >= win_48_start, collect_dttm <= win_48_end) %>%
      group_by(hospitalization_id) %>%
      summarise(any_positive_pulm_culture_48h = TRUE, .groups = "drop"),
    by = "hospitalization_id"
  ) %>%
  left_join(
    micro %>%
      inner_join(shrf_cohort %>% select(hospitalization_id, early_window_start, early_window_end), by = "hospitalization_id") %>%
      filter(collect_dttm >= early_window_start, collect_dttm <= early_window_end) %>%
      group_by(hospitalization_id) %>%
      summarise(any_positive_pulm_culture_early = TRUE, .groups = "drop"),
    by = "hospitalization_id"
  ) %>%
  mutate(
    across(starts_with("any_positive_pulm"), ~ coalesce(.x, FALSE))
  )

pm25_annual <- arrow::read_parquet(pm25_path) %>%
  transmute(
    zipcode_five_digit = normalize_zip(zip),
    admission_year = as.integer(year),
    pm25_annual = as.numeric(pm25_ug_m3)
  ) %>%
  group_by(zipcode_five_digit, admission_year) %>%
  summarise(pm25_annual = mean(pm25_annual, na.rm = TRUE), .groups = "drop")

no2_annual <- arrow::read_parquet(no2_path) %>%
  transmute(
    zipcode_five_digit = normalize_zip(zip),
    admission_year = as.integer(year),
    no2_annual = as.numeric(no2)
  ) %>%
  distinct(zipcode_five_digit, admission_year, .keep_all = TRUE)

o3_annual <- arrow::read_parquet(o3_path) %>%
  transmute(
    zipcode_five_digit = normalize_zip(zip),
    admission_year = as.integer(year),
    o3_annual = as.numeric(o3_ppb)
  ) %>%
  group_by(zipcode_five_digit, admission_year) %>%
  summarise(o3_annual = mean(o3_annual, na.rm = TRUE), .groups = "drop")

analysis_dat <- shrf_cohort %>%
  select(
    patient_id, hospitalization_id, age_at_admission, age_band, sex,
    admission_year, zipcode_five_digit
  ) %>%
  left_join(culture_flags %>% select(hospitalization_id, starts_with("any_positive_pulm")), by = "hospitalization_id") %>%
  left_join(pm25_annual, by = c("zipcode_five_digit", "admission_year")) %>%
  left_join(o3_annual, by = c("zipcode_five_digit", "admission_year")) %>%
  left_join(no2_annual, by = c("zipcode_five_digit", "admission_year"))

model_grid <- tidyr::expand_grid(
  outcome = "any_positive_pulm_culture_48h",
  exposure = c("pm25_annual", "o3_annual", "no2_annual")
)

model_results <- purrr::pmap_dfr(model_grid, ~ fit_pollution_model(analysis_dat, outcome = ..1, exposure = ..2)) %>%
  mutate(fdr_p_value = p.adjust(p_value, method = "BH"))

coverage_summary <- tibble(
  n_shrf_hospitalizations = nrow(analysis_dat),
  n_shrf_patients = n_distinct(analysis_dat$patient_id),
  n_with_zip = sum(!is.na(analysis_dat$zipcode_five_digit)),
  n_with_pm25 = sum(!is.na(analysis_dat$pm25_annual)),
  n_with_o3 = sum(!is.na(analysis_dat$o3_annual)),
  n_with_no2 = sum(!is.na(analysis_dat$no2_annual)),
  n_any_positive_pulm_culture_48h = sum(analysis_dat$any_positive_pulm_culture_48h, na.rm = TRUE)
)

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_path <- file.path(out_dir, glue("shrf_zcta_pollution_any_pulmonary_culture_models_{site_name}_{stamp}.csv"))
coverage_path <- file.path(out_dir, glue("shrf_zcta_pollution_any_pulmonary_culture_coverage_{site_name}_{stamp}.csv"))

readr::write_csv(model_results, model_path)
readr::write_csv(coverage_summary, coverage_path)

message("Coverage:")
print(coverage_summary)
message("")
message("Models:")
print(model_results, n = nrow(model_results))
message("Wrote model results: ", model_path)
message("Wrote coverage summary: ", coverage_path)
