# ================================================================================================
# SHRF ZCTA Pollution vs Specific Pulmonary Organism Models
#
# Cohort:
#   Adult ED -> ICU hospitalizations with severe hypoxemic respiratory failure:
#   IMV + PaO2/FiO2 < 300 in first 24h after first ICU admission.
#   ED intubations are eligible.
#
# Exposure:
#   ZCTA-level PM2.5 and ozone annualized from monthly data, and annual NO2.
#
# Outcomes:
#   Organism-specific positive pulmonary culture indicators.
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

MIN_ORGANISM_DETECTIONS <- as.integer(Sys.getenv("MIN_ORGANISM_DETECTIONS", unset = "25"))
WINDOWS_TO_MODEL <- str_split(Sys.getenv("CULTURE_WINDOWS", unset = "first_48h"), ",", simplify = TRUE) %>%
  as.character() %>%
  str_trim()

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
  cut(age, breaks = c(18, 40, 65, 75, Inf), right = FALSE, labels = c("18-39", "40-64", "65-74", "75+"))
}

harmonize_sex <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x %in% c("female", "f") ~ "Female",
    x %in% c("male", "m") ~ "Male",
    TRUE ~ "Other/Unknown"
  )
}

fit_pollution_model <- function(data, organism, exposure, window) {
  outcome_col <- paste(organism, window, sep = "__")
  dat <- data %>%
    mutate(
      outcome = as.integer(.data[[outcome_col]]),
      exposure_value = .data[[exposure]]
    ) %>%
    filter(!is.na(outcome), !is.na(exposure_value), !is.na(zipcode_five_digit), !is.na(admission_year))

  if (sum(dat$outcome, na.rm = TRUE) < 20) return(NULL)

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

  fit <- tryCatch(
    glmmTMB(
      outcome ~ exposure_iqr_scaled + age_band + sex + admission_year + (1 | zipcode_five_digit),
      family = binomial(),
      data = dat,
      control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
    ),
    error = function(e) {
      warning("Model failed for ", organism, " / ", exposure, " / ", window, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  coef_tab <- summary(fit)$coefficients$cond
  beta <- coef_tab["exposure_iqr_scaled", "Estimate"]
  se <- coef_tab["exposure_iqr_scaled", "Std. Error"]
  p_value <- coef_tab["exposure_iqr_scaled", "Pr(>|z|)"]

  tibble(
    organism_category = organism,
    window = window,
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
message("Minimum organism detections: ", MIN_ORGANISM_DETECTIONS)
message("Culture windows to model: ", paste(WINDOWS_TO_MODEL, collapse = ", "))

patient <- read_tbl("patient") %>% transmute(patient_id, sex_category)

hospitalization <- read_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    admission_year = year(admission_dttm),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission)),
    zipcode_five_digit = normalize_zip(zipcode_five_digit)
  ) %>%
  left_join(patient, by = "patient_id") %>%
  mutate(age_band = age_band_4(age_at_admission), sex = harmonize_sex(sex_category))

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
  summarise(first_icu_in = min(in_dttm), .groups = "drop")

pathway_flags <- adt %>%
  inner_join(icu_bounds, by = "hospitalization_id") %>%
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
  summarise(any_imv_24h = any(device_category %in% c("imv", "invasive_mechanical_ventilation"), na.rm = TRUE), .groups = "drop")

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
    positive_culture = !is.na(organism_category) & !no_growth,
    pulmonary_primary = fluid_category %in% pulmonary_primary
  ) %>%
  filter(method_category == "culture", positive_culture, pulmonary_primary)

top_organisms <- micro %>%
  distinct(hospitalization_id, organism_category) %>%
  count(organism_category, sort = TRUE, name = "n_events_any_hosp") %>%
  filter(n_events_any_hosp >= MIN_ORGANISM_DETECTIONS) %>%
  pull(organism_category)

message("Modeling ", length(top_organisms), " organisms")

make_window_wide <- function(window_name, start_col = NULL, end_col = NULL) {
  dat <- micro %>%
    filter(organism_category %in% top_organisms) %>%
    inner_join(
      shrf_cohort %>% select(hospitalization_id, all_of(c(start_col, end_col))),
      by = "hospitalization_id"
    )

  if (!is.null(start_col) && !is.null(end_col)) {
    dat <- dat %>% filter(collect_dttm >= .data[[start_col]], collect_dttm <= .data[[end_col]])
  }

  dat %>%
    distinct(hospitalization_id, organism_category) %>%
    mutate(value = 1L, organism_window = paste(organism_category, window_name, sep = "__")) %>%
    select(hospitalization_id, organism_window, value) %>%
    pivot_wider(names_from = organism_window, values_from = value, values_fill = 0L)
}

organism_wide <- shrf_cohort %>%
  select(hospitalization_id) %>%
  left_join(make_window_wide("any_hosp"), by = "hospitalization_id") %>%
  left_join(make_window_wide("first_24h", "win_24_start", "win_24_end"), by = "hospitalization_id") %>%
  left_join(make_window_wide("first_48h", "win_48_start", "win_48_end"), by = "hospitalization_id") %>%
  left_join(make_window_wide("early", "early_window_start", "early_window_end"), by = "hospitalization_id") %>%
  mutate(across(-hospitalization_id, ~ coalesce(.x, 0L)))

pm25_annual <- arrow::read_parquet(pm25_path) %>%
  transmute(zipcode_five_digit = normalize_zip(zip), admission_year = as.integer(year), pm25_annual = as.numeric(pm25_ug_m3)) %>%
  group_by(zipcode_five_digit, admission_year) %>%
  summarise(pm25_annual = mean(pm25_annual, na.rm = TRUE), .groups = "drop")

no2_annual <- arrow::read_parquet(no2_path) %>%
  transmute(zipcode_five_digit = normalize_zip(zip), admission_year = as.integer(year), no2_annual = as.numeric(no2)) %>%
  distinct(zipcode_five_digit, admission_year, .keep_all = TRUE)

o3_annual <- arrow::read_parquet(o3_path) %>%
  transmute(zipcode_five_digit = normalize_zip(zip), admission_year = as.integer(year), o3_annual = as.numeric(o3_ppb)) %>%
  group_by(zipcode_five_digit, admission_year) %>%
  summarise(o3_annual = mean(o3_annual, na.rm = TRUE), .groups = "drop")

analysis_dat <- shrf_cohort %>%
  select(patient_id, hospitalization_id, age_at_admission, age_band, sex, admission_year, zipcode_five_digit) %>%
  left_join(organism_wide, by = "hospitalization_id") %>%
  left_join(pm25_annual, by = c("zipcode_five_digit", "admission_year")) %>%
  left_join(o3_annual, by = c("zipcode_five_digit", "admission_year")) %>%
  left_join(no2_annual, by = c("zipcode_five_digit", "admission_year"))

model_grid <- tidyr::expand_grid(
  organism = top_organisms,
  window = WINDOWS_TO_MODEL,
  exposure = c("pm25_annual", "o3_annual", "no2_annual")
)

model_results <- purrr::pmap_dfr(model_grid, ~ fit_pollution_model(analysis_dat, organism = ..1, exposure = ..3, window = ..2)) %>%
  group_by(window, exposure) %>%
  mutate(fdr_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(window, exposure, fdr_p_value, desc(odds_ratio_per_iqr))

organism_counts <- micro %>%
  distinct(hospitalization_id, organism_category) %>%
  count(organism_category, sort = TRUE, name = "n_shrf_hospitalizations_any_hosp")

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_path <- file.path(out_dir, glue("shrf_zcta_pollution_organism_models_{site_name}_{stamp}.csv"))
counts_path <- file.path(out_dir, glue("shrf_zcta_pollution_organism_model_counts_{site_name}_{stamp}.csv"))

readr::write_csv(model_results, model_path)
readr::write_csv(organism_counts, counts_path)

message("Top positive associations, FDR < 0.10:")
model_results %>%
  filter(odds_ratio_per_iqr > 1, fdr_p_value < 0.10) %>%
  select(organism_category, window, exposure, n_events, n_hospitalizations, n_zctas,
         odds_ratio_per_iqr, ci_low, ci_high, p_value, fdr_p_value) %>%
  print(n = 100)

message("Wrote model results: ", model_path)
message("Wrote organism counts: ", counts_path)
