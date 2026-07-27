# ================================================================================================
# ED-to-ICU Early Respiratory Culture ZCTA Pollution vs Organism Models
# Cohort:
#   Adult ED -> ICU hospitalizations with, in the first 48h after hospital admission:
#     1) positive pulmonary respiratory culture
#     2) oxygen >= 2 L/min or advanced oxygen/ventilatory support
#     3) administered antibacterial antibiotic
#
# Outcomes:
#   Organism-specific positive pulmonary culture indicators in the first 48h after admission.
# ================================================================================================

suppressPackageStartupMessages({
  library(arrow)
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

fit_pollution_model <- function(data, organism, exposure) {
  dat <- data %>%
    mutate(
      outcome = as.integer(.data[[organism]]),
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
      warning("Model failed for ", organism, " / ", exposure, ": ", conditionMessage(e))
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

pulmonary_primary <- c("respiratory_tract", "respiratory_tract_lower")
advanced_o2_devices <- c("imv", "invasive_mechanical_ventilation", "nippv", "niv", "cpap", "high flow nc", "high_flow_nasal_cannula", "hfnc", "face mask", "trach collar")
antibacterial_categories <- c(
  "amikacin", "amoxicillin", "amoxicillin_clavulanate", "ampicillin", "ampicillin_sulbactam",
  "azithromycin", "aztreonam", "cefadroxil", "cefazolin", "cefdinir", "cefepime", "cefixime",
  "cefotaxime", "cefoxitin", "cefpodoxime", "ceftaroline", "ceftazidime", "ceftazidime_avibactam",
  "ceftriaxone", "cefuroxime", "cephalexin", "ciprofloxacin", "clarithromycin", "clindamycin",
  "daptomycin", "dicloxacillin", "doxycycline", "ertapenem", "erythromycin", "fidaxomicin",
  "fosfomycin", "gentamicin", "imipenem", "imipenem_relebactam", "levofloxacin", "linezolid",
  "meropenem", "metronidazole", "minocycline", "moxifloxacin", "nitrofurantoin", "oxacillin",
  "penicillin", "piperacillin_tazobactam", "quinupristin_dalfopristin", "rifampin", "streptomycin",
  "sulfadiazine", "sulbactam_durlobactam", "tedizolid", "tetracycline", "tigecycline", "tobramycin",
  "trimethoprim", "trimethoprim_sulfamethoxazole", "vancomycin"
)

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
    location_category = str_to_lower(str_trim(as.character(location_category)))
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
    window_start = admission_dttm,
    window_end = admission_dttm + hours(48)
  ) %>%
  filter(adult, ed_to_icu, !is.na(window_start))

base_ids <- base$hospitalization_id

micro <- read_tbl("microbiology_culture") %>%
  filter(hospitalization_id %in% base_ids) %>%
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
  filter(method_category == "culture", pulmonary_primary) %>%
  inner_join(base %>% select(hospitalization_id, window_start, window_end), by = "hospitalization_id") %>%
  filter(!is.na(collect_dttm), collect_dttm >= window_start, collect_dttm <= window_end)

positive_culture_flags <- micro %>%
  group_by(hospitalization_id) %>%
  summarise(any_positive_resp_culture_48h = any(positive_culture, na.rm = TRUE), .groups = "drop")

resp_support_flags <- read_tbl("respiratory_support") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    device_category = str_to_lower(str_trim(as.character(device_category))),
    lpm_set = suppressWarnings(as.numeric(lpm_set))
  ) %>%
  inner_join(base %>% select(hospitalization_id, window_start, window_end), by = "hospitalization_id") %>%
  filter(!is.na(recorded_dttm), recorded_dttm >= window_start, recorded_dttm <= window_end) %>%
  group_by(hospitalization_id) %>%
  summarise(any_o2_ge2l_or_advanced_support_48h = any(lpm_set >= 2 | device_category %in% advanced_o2_devices, na.rm = TRUE), .groups = "drop")

read_meds <- function(table_name) {
  tbl <- read_tbl(table_name, required = FALSE)
  if (is.null(tbl)) {
    return(tibble(hospitalization_id = character(), admin_dttm = as.POSIXct(character()), med_category = character(), mar_action_group = character()))
  }
  tbl %>%
    transmute(
      hospitalization_id,
      admin_dttm = safe_ts(admin_dttm),
      med_category = str_to_lower(str_trim(as.character(med_category))),
      mar_action_group = str_to_lower(str_trim(as.character(mar_action_group)))
    )
}

antibiotic_flags <- bind_rows(
  read_meds("medication_admin_intermittent"),
  read_meds("medication_admin_continuous")
) %>%
  filter(hospitalization_id %in% base_ids) %>%
  inner_join(base %>% select(hospitalization_id, window_start, window_end), by = "hospitalization_id") %>%
  filter(!is.na(admin_dttm), admin_dttm >= window_start, admin_dttm <= window_end) %>%
  mutate(
    administered = is.na(mar_action_group) | mar_action_group == "administered",
    antibacterial_antibiotic = med_category %in% antibacterial_categories
  ) %>%
  group_by(hospitalization_id) %>%
  summarise(any_antibacterial_antibiotic_48h = any(antibacterial_antibiotic & administered, na.rm = TRUE), .groups = "drop")

analysis_cohort <- base %>%
  left_join(positive_culture_flags, by = "hospitalization_id") %>%
  left_join(resp_support_flags, by = "hospitalization_id") %>%
  left_join(antibiotic_flags, by = "hospitalization_id") %>%
  mutate(
    any_positive_resp_culture_48h = coalesce(any_positive_resp_culture_48h, FALSE),
    any_o2_ge2l_or_advanced_support_48h = coalesce(any_o2_ge2l_or_advanced_support_48h, FALSE),
    any_antibacterial_antibiotic_48h = coalesce(any_antibacterial_antibiotic_48h, FALSE)
  ) %>%
  filter(any_positive_resp_culture_48h, any_o2_ge2l_or_advanced_support_48h, any_antibacterial_antibiotic_48h)

top_organisms <- micro %>%
  filter(positive_culture, hospitalization_id %in% analysis_cohort$hospitalization_id) %>%
  distinct(hospitalization_id, organism_category) %>%
  count(organism_category, sort = TRUE, name = "n_events") %>%
  filter(n_events >= MIN_ORGANISM_DETECTIONS) %>%
  pull(organism_category)

message("Analysis cohort hospitalizations: ", nrow(analysis_cohort))
message("Modeling ", length(top_organisms), " organisms")

organism_wide <- micro %>%
  filter(positive_culture, hospitalization_id %in% analysis_cohort$hospitalization_id, organism_category %in% top_organisms) %>%
  distinct(hospitalization_id, organism_category) %>%
  mutate(value = 1L) %>%
  pivot_wider(names_from = organism_category, values_from = value, values_fill = 0L)

pm25_annual <- arrow::read_parquet(pm25_path) %>%
  transmute(zipcode_five_digit = normalize_zip(zip), admission_year = as.integer(year), pm25_annual = as.numeric(pm25_ug_m3)) %>%
  group_by(zipcode_five_digit, admission_year) %>%
  summarise(pm25_annual = mean(pm25_annual, na.rm = TRUE), .groups = "drop")

o3_annual <- arrow::read_parquet(o3_path) %>%
  transmute(zipcode_five_digit = normalize_zip(zip), admission_year = as.integer(year), o3_annual = as.numeric(o3_ppb)) %>%
  group_by(zipcode_five_digit, admission_year) %>%
  summarise(o3_annual = mean(o3_annual, na.rm = TRUE), .groups = "drop")

no2_annual <- arrow::read_parquet(no2_path) %>%
  transmute(zipcode_five_digit = normalize_zip(zip), admission_year = as.integer(year), no2_annual = as.numeric(no2)) %>%
  distinct(zipcode_five_digit, admission_year, .keep_all = TRUE)

analysis_dat <- analysis_cohort %>%
  select(patient_id, hospitalization_id, age_at_admission, age_band, sex, admission_year, zipcode_five_digit) %>%
  left_join(organism_wide, by = "hospitalization_id") %>%
  mutate(across(all_of(top_organisms), ~ coalesce(.x, 0L))) %>%
  left_join(pm25_annual, by = c("zipcode_five_digit", "admission_year")) %>%
  left_join(o3_annual, by = c("zipcode_five_digit", "admission_year")) %>%
  left_join(no2_annual, by = c("zipcode_five_digit", "admission_year"))

model_grid <- tidyr::expand_grid(
  organism = top_organisms,
  exposure = c("pm25_annual", "o3_annual", "no2_annual")
)

model_results <- purrr::pmap_dfr(model_grid, ~ fit_pollution_model(analysis_dat, organism = ..1, exposure = ..2)) %>%
  group_by(exposure) %>%
  mutate(fdr_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(exposure, fdr_p_value, desc(odds_ratio_per_iqr))

organism_counts <- micro %>%
  filter(positive_culture, hospitalization_id %in% analysis_cohort$hospitalization_id) %>%
  distinct(hospitalization_id, organism_category) %>%
  count(organism_category, sort = TRUE, name = "n_cohort_hospitalizations")

coverage_summary <- tibble(
  n_cohort_hospitalizations = nrow(analysis_dat),
  n_cohort_patients = n_distinct(analysis_dat$patient_id),
  n_with_zip = sum(!is.na(analysis_dat$zipcode_five_digit)),
  n_with_pm25 = sum(!is.na(analysis_dat$pm25_annual)),
  n_with_o3 = sum(!is.na(analysis_dat$o3_annual)),
  n_with_no2 = sum(!is.na(analysis_dat$no2_annual)),
  n_modeled_organisms = length(top_organisms)
)

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_path <- file.path(out_dir, glue("ed_icu_early_resp_culture_pollution_organism_models_{site_name}_{stamp}.csv"))
counts_path <- file.path(out_dir, glue("ed_icu_early_resp_culture_organism_counts_{site_name}_{stamp}.csv"))
coverage_path <- file.path(out_dir, glue("ed_icu_early_resp_culture_pollution_coverage_{site_name}_{stamp}.csv"))

readr::write_csv(model_results, model_path)
readr::write_csv(organism_counts, counts_path)
readr::write_csv(coverage_summary, coverage_path)

message("Coverage:")
print(coverage_summary)
message("")
message("Top positive associations, FDR < 0.10:")
model_results %>%
  filter(odds_ratio_per_iqr > 1, fdr_p_value < 0.10) %>%
  select(organism_category, exposure, n_events, n_hospitalizations, n_zctas,
         odds_ratio_per_iqr, ci_low, ci_high, p_value, fdr_p_value) %>%
  print(n = 100)
message("Wrote model results: ", model_path)
message("Wrote organism counts: ", counts_path)
message("Wrote coverage summary: ", coverage_path)
