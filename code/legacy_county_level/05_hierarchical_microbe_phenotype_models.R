# ================================================================================================
# CLIF Pollution-Microbiome | Patient-Level Hierarchical Phenotype Models
# Purpose:
#   Fit mixed-effects logistic models testing whether PM2.5 or NO2 correspond to higher odds of
#   detecting specific respiratory organisms, allowing effects to differ by pneumonia, sepsis, and
#   severe respiratory support phenotypes.
#
# Model family:
#   organism_detected ~ exposure_iqr_scaled * phenotype + age_band + sex + year + (1 | county_fips)
#
# Denominator:
#   Adult ICU hospitalizations with at least one primary respiratory culture in the early ICU window.
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

source("utils/config.R")

`%||%` <- function(x, y) if (!is.null(x)) x else y

site_name <- config$site_name
tables_path <- config$tables_path
exposome_path <- config$exposome_path %||% file.path(config$repo %||% getwd(), "exposome")
analysis_label <- Sys.getenv("ANALYSIS_LABEL", unset = "all_counties")
exclude_county_fips_raw <- Sys.getenv("EXCLUDE_COUNTY_FIPS", unset = "")

START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
END_DATE <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")
MIN_ICU_LOS_H <- 24
MICROBE_WINDOW_PRE_H <- 48
MICROBE_WINDOW_POST_H <- 72
MIN_ORGANISM_DETECTIONS <- as.integer(Sys.getenv("MIN_ORGANISM_DETECTIONS", unset = "75"))
MAX_ORGANISMS <- as.integer(Sys.getenv("MAX_ORGANISMS", unset = "20"))
resp_fluid_mode <- Sys.getenv("RESP_FLUID_MODE", unset = "primary")
RESP_FLUID_PRIMARY <- if (identical(resp_fluid_mode, "lower_only")) {
  c("respiratory_tract_lower")
} else {
  c("respiratory_tract", "respiratory_tract_lower")
}

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
      if (!requireNamespace("fst", quietly = TRUE)) stop("Package 'fst' is required to read fst files.")
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

exclude_county_fips <- if (nzchar(exclude_county_fips_raw)) {
  str_split(exclude_county_fips_raw, ",", simplify = TRUE) %>%
    as.character() %>%
    str_trim() %>%
    normalize_county_fips()
} else {
  character()
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

extract_effect <- function(fit, exposure_term = "exposure_iqr_scaled", interaction_term = "exposure_iqr_scaled:phenotype_valueTRUE") {
  coef_tab <- summary(fit)$coefficients$cond
  vc <- vcov(fit)$cond

  if (!exposure_term %in% rownames(coef_tab)) return(NULL)

  beta_absent <- coef_tab[exposure_term, "Estimate"]
  se_absent <- coef_tab[exposure_term, "Std. Error"]
  p_absent <- coef_tab[exposure_term, "Pr(>|z|)"]

  has_interaction <- interaction_term %in% rownames(coef_tab)
  if (has_interaction) {
    beta_present <- beta_absent + coef_tab[interaction_term, "Estimate"]
    se_present <- sqrt(
      vc[exposure_term, exposure_term] +
        vc[interaction_term, interaction_term] +
        2 * vc[exposure_term, interaction_term]
    )
    z_present <- beta_present / se_present
    p_present <- 2 * pnorm(abs(z_present), lower.tail = FALSE)
    interaction_or <- exp(coef_tab[interaction_term, "Estimate"])
    interaction_p <- coef_tab[interaction_term, "Pr(>|z|)"]
  } else {
    beta_present <- NA_real_
    se_present <- NA_real_
    p_present <- NA_real_
    interaction_or <- NA_real_
    interaction_p <- NA_real_
  }

  bind_rows(
    tibble(
      phenotype_level = "absent",
      log_or = beta_absent,
      se = se_absent,
      odds_ratio_per_iqr = exp(beta_absent),
      ci_low = exp(beta_absent - 1.96 * se_absent),
      ci_high = exp(beta_absent + 1.96 * se_absent),
      p_value = p_absent,
      interaction_or = interaction_or,
      interaction_p_value = interaction_p
    ),
    tibble(
      phenotype_level = "present",
      log_or = beta_present,
      se = se_present,
      odds_ratio_per_iqr = exp(beta_present),
      ci_low = exp(beta_present - 1.96 * se_present),
      ci_high = exp(beta_present + 1.96 * se_present),
      p_value = p_present,
      interaction_or = interaction_or,
      interaction_p_value = interaction_p
    )
  )
}

fit_one <- function(dat, organism, exposure, phenotype) {
  exposure_iqr <- IQR(dat[[exposure]], na.rm = TRUE)
  if (!is.finite(exposure_iqr) || exposure_iqr <= 0) return(NULL)

  model_dat <- dat %>%
    transmute(
      outcome = as.integer(.data[[organism]]),
      exposure_iqr_scaled = .data[[exposure]] / exposure_iqr,
      phenotype_value = as.logical(.data[[phenotype]]),
      county_fips = factor(county_fips),
      year = factor(year),
      age_band = fct_na_value_to_level(factor(age_band), level = "Unknown"),
      sex = factor(sex)
    ) %>%
    filter(!is.na(outcome), !is.na(exposure_iqr_scaled), !is.na(phenotype_value), !is.na(county_fips))

  if (
    nrow(model_dat) < 200 ||
      sum(model_dat$outcome == 1, na.rm = TRUE) < 20 ||
      n_distinct(model_dat$county_fips) < 5 ||
      length(unique(model_dat$phenotype_value)) < 2
  ) {
    return(NULL)
  }

  fit <- tryCatch(
    glmmTMB(
      outcome ~ exposure_iqr_scaled * phenotype_value + age_band + sex + year + (1 | county_fips),
      family = binomial(),
      data = model_dat,
      control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
    ),
    error = function(e) {
      warning("Model failed for ", organism, " / ", exposure, " / ", phenotype, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  effects <- extract_effect(fit)
  if (is.null(effects)) return(NULL)

  effects %>%
    mutate(
      organism_category = organism,
      exposure = exposure,
      phenotype = phenotype,
      exposure_iqr = exposure_iqr,
      n_hospitalizations = nrow(model_dat),
      n_detected = sum(model_dat$outcome == 1, na.rm = TRUE),
      n_counties = n_distinct(model_dat$county_fips),
      n_phenotype_present = sum(model_dat$phenotype_value, na.rm = TRUE),
      .before = 1
    )
}

message(glue("Building patient-level cohort for {site_name} ({analysis_label})"))
message("Respiratory specimen mode: ", resp_fluid_mode, " [", paste(RESP_FLUID_PRIMARY, collapse = ", "), "]")
message("Minimum organism detections: ", MIN_ORGANISM_DETECTIONS)
if (length(exclude_county_fips) > 0) {
  message("Excluding county FIPS: ", paste(exclude_county_fips, collapse = ", "))
}

patient <- get_tbl("patient") %>%
  transmute(patient_id, sex_category)

hospitalization <- get_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission)),
    county_fips = normalize_county_fips(county_code)
  )

adt <- get_tbl("adt") %>%
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
      first_icu_in >= START_DATE & first_icu_in <= END_DATE &
      !(county_fips %in% exclude_county_fips)
  )

base_ids <- base %>% filter(include_base) %>% pull(hospitalization_id)

micro_win <- get_tbl("microbiology_culture") %>%
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
    positive_culture = !is.na(organism_category) & !no_growth
  ) %>%
  inner_join(
    base %>%
      filter(include_base) %>%
      transmute(
        hospitalization_id,
        culture_window_start = first_icu_in - hours(MICROBE_WINDOW_PRE_H),
        culture_window_end = first_icu_in + hours(MICROBE_WINDOW_POST_H)
      ),
    by = "hospitalization_id"
  ) %>%
  filter(
    !is.na(collect_dttm),
    collect_dttm >= culture_window_start,
    collect_dttm <= culture_window_end,
    fluid_category %in% RESP_FLUID_PRIMARY,
    method_category == "culture"
  )

micro_features <- micro_win %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_resp_culture = TRUE,
    any_positive_resp_culture = any(positive_culture, na.rm = TRUE),
    .groups = "drop"
  )

top_organisms <- micro_win %>%
  filter(positive_culture) %>%
  distinct(hospitalization_id, organism_category) %>%
  count(organism_category, name = "n_detected", sort = TRUE) %>%
  filter(n_detected >= MIN_ORGANISM_DETECTIONS) %>%
  slice_head(n = MAX_ORGANISMS) %>%
  pull(organism_category)

message(glue("Modeling {length(top_organisms)} organisms with >= {MIN_ORGANISM_DETECTIONS} detections"))

organism_wide <- micro_win %>%
  filter(positive_culture, organism_category %in% top_organisms) %>%
  distinct(hospitalization_id, organism_category) %>%
  mutate(detected = 1L) %>%
  tidyr::pivot_wider(
    names_from = organism_category,
    values_from = detected,
    values_fill = 0L,
    names_repair = "minimal"
  )

resp_support <- get_tbl("respiratory_support") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    device_category = str_to_lower(str_trim(as.character(device_category)))
  )

resp_features <- resp_support %>%
  inner_join(base %>% filter(include_base) %>% select(hospitalization_id, first_icu_in), by = "hospitalization_id") %>%
  filter(recorded_dttm >= first_icu_in - hours(24), recorded_dttm <= first_icu_in + hours(72)) %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_hfnc = any(device_category %in% c("hfnc", "high flow nc", "high_flow_nasal_cannula"), na.rm = TRUE),
    any_niv = any(device_category %in% c("niv", "nippv", "bipap", "cpap"), na.rm = TRUE),
    any_imv = any(device_category %in% c("imv", "invasive_mechanical_ventilation"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(severe_arf_proxy = any_imv | any_niv | any_hfnc)

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
  diagnosis_flags <- tibble(hospitalization_id = base_ids, pneumonia_dx = NA, sepsis_dx = NA)
}

exposure <- load_pollutant("pm25") %>%
  full_join(load_pollutant("no2"), by = c("county_fips", "year"))

model_dat <- base %>%
  filter(include_base) %>%
  left_join(micro_features, by = "hospitalization_id") %>%
  filter(coalesce(any_resp_culture, FALSE)) %>%
  left_join(resp_features, by = "hospitalization_id") %>%
  left_join(diagnosis_flags, by = "hospitalization_id") %>%
  left_join(exposure, by = c("county_fips", "year")) %>%
  left_join(organism_wide, by = "hospitalization_id") %>%
  mutate(
    across(all_of(top_organisms), ~ coalesce(.x, 0L)),
    severe_arf_proxy = coalesce(severe_arf_proxy, FALSE),
    pneumonia_dx = coalesce(pneumonia_dx, FALSE),
    sepsis_dx = coalesce(sepsis_dx, FALSE)
  )

model_grid <- tidyr::expand_grid(
  organism = top_organisms,
  exposure = c("pm25_mean", "no2_mean"),
  phenotype = c("pneumonia_dx", "sepsis_dx", "severe_arf_proxy")
)

message(glue("Fitting {nrow(model_grid)} mixed models"))

model_results <- purrr::pmap_dfr(
  model_grid,
  ~ fit_one(model_dat, organism = ..1, exposure = ..2, phenotype = ..3)
) %>%
  group_by(exposure, phenotype, phenotype_level) %>%
  mutate(fdr_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(exposure, phenotype, phenotype_level, fdr_p_value, desc(odds_ratio_per_iqr))

out_path <- file.path(out_dir, glue("hierarchical_pollution_microbe_phenotype_models_{site_name}_{analysis_label}_{stamp}.csv"))
readr::write_csv(model_results, out_path)

message("Wrote hierarchical model results: ", out_path)
message("")
message("Phenotype-present increased risk signals, FDR < 0.10:")
model_results %>%
  filter(phenotype_level == "present", odds_ratio_per_iqr > 1, fdr_p_value < 0.10) %>%
  arrange(fdr_p_value, desc(odds_ratio_per_iqr)) %>%
  select(organism_category, exposure, phenotype, n_detected, n_phenotype_present,
         odds_ratio_per_iqr, ci_low, ci_high, p_value, fdr_p_value, interaction_or, interaction_p_value) %>%
  print(n = 100)
