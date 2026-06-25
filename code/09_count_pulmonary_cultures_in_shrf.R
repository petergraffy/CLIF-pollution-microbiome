# ================================================================================================
# Count positive pulmonary cultures among severe hypoxemic respiratory failure hospitalizations
#
# Severe hypoxemic RF definition:
#   Adult ED -> ICU hospitalizations with IMV and PaO2/FiO2 < 300 in first 24h
#   after first ICU admission. ED intubations are eligible.
#
# Pulmonary culture definition:
#   microbiology_culture rows with method_category == "culture", positive organism, and
#   fluid_category in pulmonary/respiratory categories.
# ================================================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(glue)
  library(janitor)
  library(jsonlite)
  library(lubridate)
  library(readr)
  library(stringr)
})

config_path <- file.path("config", "config.json")
config <- if (file.exists(config_path)) jsonlite::fromJSON(config_path) else list()

candidate_paths <- c(
  Sys.getenv("CLIF_TABLES_PATH", unset = NA_character_),
  config$tables_path,
  "/Users/saborpete/Desktop/Peter/Postdoc/CLIF v2.1/2.1.0",
  "/Users/saborpete/Desktop/Peter/Postdoc/CLIF v2.1",
  "/Users/saborpete/Library/CloudStorage/Box-Box/03-CLIF-2.1/2.1.0"
)

tables_path <- candidate_paths[!is.na(candidate_paths) & dir.exists(candidate_paths)][1]
if (is.na(tables_path)) stop("Could not locate a CLIF table directory.")

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

read_tbl <- function(name) {
  path <- file.path(tables_path, paste0("clif_", name, ".parquet"))
  if (!file.exists(path)) stop("Missing table: ", path)
  arrow::read_parquet(path) %>% janitor::clean_names()
}

normalize_fio2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  case_when(
    x > 1 & x <= 100 ~ x / 100,
    x >= 0.21 & x <= 1 ~ x,
    TRUE ~ NA_real_
  )
}

message("Using CLIF tables: ", tables_path)

hospitalization <- read_tbl("hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission))
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
  select(
    patient_id, hospitalization_id, admission_dttm, discharge_dttm,
    first_icu_in, win_24_start, win_24_end, win_48_start, win_48_end,
    early_window_start, early_window_end
  ) %>%
  left_join(imv_flags, by = "hospitalization_id") %>%
  left_join(pf_flags, by = "hospitalization_id") %>%
  mutate(
    any_imv_24h = coalesce(any_imv_24h, FALSE),
    any_pf_lt_300_24h = coalesce(any_pf_lt_300_24h, FALSE),
    severe_hypoxemic_rf_24h = any_imv_24h & any_pf_lt_300_24h
  ) %>%
  filter(severe_hypoxemic_rf_24h)

pulmonary_primary <- c("respiratory_tract", "respiratory_tract_lower")
pulmonary_with_pleural <- c(pulmonary_primary, "pleural_cavity_fluid")

micro <- read_tbl("microbiology_culture") %>%
  filter(hospitalization_id %in% shrf_cohort$hospitalization_id) %>%
  transmute(
    patient_id,
    hospitalization_id,
    collect_dttm = safe_ts(collect_dttm),
    result_dttm = safe_ts(result_dttm),
    fluid_category = str_to_lower(str_trim(as.character(fluid_category))),
    method_category = str_to_lower(str_trim(as.character(method_category))),
    organism_category = str_to_lower(str_trim(as.character(organism_category))),
    organism_group = str_to_lower(str_trim(as.character(organism_group)))
  ) %>%
  mutate(
    organism_group = coalesce(na_if(organism_group, ""), organism_category),
    no_growth = organism_group %in% c("no_growth", "no growth"),
    positive_culture = !is.na(organism_group) & !no_growth,
    pulmonary_primary = fluid_category %in% pulmonary_primary,
    pulmonary_with_pleural = fluid_category %in% pulmonary_with_pleural
  ) %>%
  filter(method_category == "culture", positive_culture)

count_window <- function(data, cohort, fluid_col, window_name, start_col = NULL, end_col = NULL) {
  dat <- data %>%
    filter(.data[[fluid_col]]) %>%
    inner_join(
      cohort %>% select(patient_id, hospitalization_id, all_of(c(start_col, end_col))),
      by = c("patient_id", "hospitalization_id")
    )

  if (!is.null(start_col) && !is.null(end_col)) {
    dat <- dat %>% filter(!is.na(collect_dttm), collect_dttm >= .data[[start_col]], collect_dttm <= .data[[end_col]])
  }

  tibble(
    definition = paste(fluid_col, window_name, sep = "__"),
    n_shrf_hospitalizations = n_distinct(cohort$hospitalization_id),
    n_shrf_patients = n_distinct(cohort$patient_id),
    n_positive_pulmonary_culture_hospitalizations = n_distinct(dat$hospitalization_id),
    n_positive_pulmonary_culture_patients = n_distinct(dat$patient_id)
  )
}

summary_tbl <- bind_rows(
  count_window(micro, shrf_cohort, "pulmonary_primary", "any_time_during_hospitalization"),
  count_window(micro, shrf_cohort, "pulmonary_primary", "first_icu_24h", "win_24_start", "win_24_end"),
  count_window(micro, shrf_cohort, "pulmonary_primary", "first_icu_48h", "win_48_start", "win_48_end"),
  count_window(micro, shrf_cohort, "pulmonary_primary", "early_icu_minus48_plus72h", "early_window_start", "early_window_end"),
  count_window(micro, shrf_cohort, "pulmonary_with_pleural", "any_time_during_hospitalization"),
  count_window(micro, shrf_cohort, "pulmonary_with_pleural", "first_icu_24h", "win_24_start", "win_24_end"),
  count_window(micro, shrf_cohort, "pulmonary_with_pleural", "first_icu_48h", "win_48_start", "win_48_end"),
  count_window(micro, shrf_cohort, "pulmonary_with_pleural", "early_icu_minus48_plus72h", "early_window_start", "early_window_end")
) %>%
  mutate(
    pct_hospitalizations = 100 * n_positive_pulmonary_culture_hospitalizations / n_shrf_hospitalizations,
    pct_patients = 100 * n_positive_pulmonary_culture_patients / n_shrf_patients
  )

organism_summary <- micro %>%
  filter(pulmonary_primary) %>%
  distinct(hospitalization_id, organism_category) %>%
  count(organism_category, sort = TRUE, name = "n_shrf_hospitalizations") %>%
  slice_head(n = 25)

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
summary_path <- file.path(out_dir, glue("shrf_positive_pulmonary_culture_count_UCMC_{stamp}.csv"))
organism_path <- file.path(out_dir, glue("shrf_positive_pulmonary_culture_top_organisms_UCMC_{stamp}.csv"))

readr::write_csv(summary_tbl, summary_path)
readr::write_csv(organism_summary, organism_path)

print(summary_tbl)
message("Wrote count summary: ", summary_path)
message("Wrote top organism summary: ", organism_path)
