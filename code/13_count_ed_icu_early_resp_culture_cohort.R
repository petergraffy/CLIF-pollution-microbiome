# ================================================================================================
# Count ED-to-ICU Early Respiratory Culture Cohort
# Definition:
#   Adult ED -> ICU hospitalizations with respiratory cultures in the first 48 hours after hospital
#   admission, then sequentially filtered to:
#     1) positive respiratory culture
#     2) minimum 2 L/min oxygen or advanced oxygen/ventilatory support
#     3) antibacterial antibiotic administration
#
# Window:
#   Hospital admission through admission + 48 hours.
# ================================================================================================

suppressPackageStartupMessages({
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

message("Using CLIF tables: ", tables_path)

pulmonary_primary <- c("respiratory_tract", "respiratory_tract_lower")
advanced_o2_devices <- c(
  "imv",
  "invasive_mechanical_ventilation",
  "nippv",
  "niv",
  "cpap",
  "high flow nc",
  "high_flow_nasal_cannula",
  "hfnc",
  "face mask",
  "trach collar"
)

antibacterial_categories <- c(
  "amikacin",
  "amoxicillin",
  "amoxicillin_clavulanate",
  "ampicillin",
  "ampicillin_sulbactam",
  "azithromycin",
  "aztreonam",
  "cefadroxil",
  "cefazolin",
  "cefdinir",
  "cefepime",
  "cefixime",
  "cefotaxime",
  "cefoxitin",
  "cefpodoxime",
  "ceftaroline",
  "ceftazidime",
  "ceftazidime_avibactam",
  "ceftriaxone",
  "cefuroxime",
  "cephalexin",
  "ciprofloxacin",
  "clarithromycin",
  "clindamycin",
  "daptomycin",
  "dicloxacillin",
  "doxycycline",
  "ertapenem",
  "erythromycin",
  "fidaxomicin",
  "fosfomycin",
  "gentamicin",
  "imipenem",
  "imipenem_relebactam",
  "levofloxacin",
  "linezolid",
  "meropenem",
  "metronidazole",
  "minocycline",
  "moxifloxacin",
  "nitrofurantoin",
  "oxacillin",
  "penicillin",
  "piperacillin_tazobactam",
  "quinupristin_dalfopristin",
  "rifampin",
  "streptomycin",
  "sulfadiazine",
  "sulbactam_durlobactam",
  "tedizolid",
  "tetracycline",
  "tigecycline",
  "tobramycin",
  "trimethoprim",
  "trimethoprim_sulfamethoxazole",
  "vancomycin"
)

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
  filter(method_category == "culture", pulmonary_primary) %>%
  inner_join(
    base %>% select(hospitalization_id, window_start, window_end),
    by = "hospitalization_id"
  ) %>%
  filter(!is.na(collect_dttm), collect_dttm >= window_start, collect_dttm <= window_end)

culture_flags <- micro %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_resp_culture_48h = TRUE,
    any_positive_resp_culture_48h = any(positive_culture, na.rm = TRUE),
    n_resp_culture_rows_48h = n(),
    n_positive_resp_culture_rows_48h = sum(positive_culture, na.rm = TRUE),
    .groups = "drop"
  )

resp_support <- read_tbl("respiratory_support") %>%
  filter(hospitalization_id %in% base_ids) %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    device_category = str_to_lower(str_trim(as.character(device_category))),
    lpm_set = suppressWarnings(as.numeric(lpm_set))
  ) %>%
  inner_join(
    base %>% select(hospitalization_id, window_start, window_end),
    by = "hospitalization_id"
  ) %>%
  filter(!is.na(recorded_dttm), recorded_dttm >= window_start, recorded_dttm <= window_end)

oxygen_flags <- resp_support %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_o2_ge2l_or_advanced_support_48h = any(lpm_set >= 2 | device_category %in% advanced_o2_devices, na.rm = TRUE),
    max_lpm_set_48h = suppressWarnings(max(lpm_set, na.rm = TRUE)),
    any_advanced_o2_support_48h = any(device_category %in% advanced_o2_devices, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(max_lpm_set_48h = if_else(is.infinite(max_lpm_set_48h), NA_real_, max_lpm_set_48h))

read_meds <- function(table_name) {
  tbl <- read_tbl(table_name, required = FALSE)
  if (is.null(tbl)) {
    return(tibble(
      hospitalization_id = character(),
      admin_dttm = as.POSIXct(character()),
      med_category = character(),
      med_group = character(),
      mar_action_group = character()
    ))
  }

  tbl %>%
    transmute(
      hospitalization_id,
      admin_dttm = safe_ts(admin_dttm),
      med_category = str_to_lower(str_trim(as.character(med_category))),
      med_group = str_to_lower(str_trim(as.character(med_group))),
      mar_action_group = str_to_lower(str_trim(as.character(mar_action_group)))
    )
}

meds <- bind_rows(
  read_meds("medication_admin_intermittent"),
  read_meds("medication_admin_continuous")
) %>%
  filter(hospitalization_id %in% base_ids) %>%
  inner_join(
    base %>% select(hospitalization_id, window_start, window_end),
    by = "hospitalization_id"
  ) %>%
  filter(!is.na(admin_dttm), admin_dttm >= window_start, admin_dttm <= window_end) %>%
  mutate(
    administered = is.na(mar_action_group) | mar_action_group == "administered",
    antibacterial_antibiotic = med_category %in% antibacterial_categories
  )

antibiotic_flags <- meds %>%
  group_by(hospitalization_id) %>%
  summarise(
    any_antibacterial_antibiotic_48h = any(antibacterial_antibiotic & administered, na.rm = TRUE),
    any_cms_sepsis_qualifying_antimicrobial_48h = any(med_group == "cms_sepsis_qualifying_antibiotics" & administered, na.rm = TRUE),
    .groups = "drop"
  )

cohort_flags <- base %>%
  select(patient_id, hospitalization_id) %>%
  left_join(culture_flags, by = "hospitalization_id") %>%
  left_join(oxygen_flags, by = "hospitalization_id") %>%
  left_join(antibiotic_flags, by = "hospitalization_id") %>%
  mutate(
    any_resp_culture_48h = coalesce(any_resp_culture_48h, FALSE),
    any_positive_resp_culture_48h = coalesce(any_positive_resp_culture_48h, FALSE),
    any_o2_ge2l_or_advanced_support_48h = coalesce(any_o2_ge2l_or_advanced_support_48h, FALSE),
    any_advanced_o2_support_48h = coalesce(any_advanced_o2_support_48h, FALSE),
    any_antibacterial_antibiotic_48h = coalesce(any_antibacterial_antibiotic_48h, FALSE),
    any_cms_sepsis_qualifying_antimicrobial_48h = coalesce(any_cms_sepsis_qualifying_antimicrobial_48h, FALSE),
    final_cohort = any_resp_culture_48h &
      any_positive_resp_culture_48h &
      any_o2_ge2l_or_advanced_support_48h &
      any_antibacterial_antibiotic_48h
  )

step_summary <- tibble(
  step = c(
    "adult_ed_to_icu",
    "resp_culture_first_48h_admission",
    "positive_resp_culture_first_48h_admission",
    "positive_resp_culture_plus_o2_ge2l_or_advanced_support_first_48h",
    "positive_resp_culture_plus_o2_plus_antibacterial_antibiotic_first_48h"
  ),
  n_hospitalizations = c(
    nrow(base),
    sum(cohort_flags$any_resp_culture_48h, na.rm = TRUE),
    sum(cohort_flags$any_resp_culture_48h & cohort_flags$any_positive_resp_culture_48h, na.rm = TRUE),
    sum(cohort_flags$any_resp_culture_48h & cohort_flags$any_positive_resp_culture_48h &
          cohort_flags$any_o2_ge2l_or_advanced_support_48h, na.rm = TRUE),
    sum(cohort_flags$final_cohort, na.rm = TRUE)
  ),
  n_patients = c(
    n_distinct(base$patient_id),
    n_distinct(cohort_flags$patient_id[cohort_flags$any_resp_culture_48h]),
    n_distinct(cohort_flags$patient_id[cohort_flags$any_resp_culture_48h & cohort_flags$any_positive_resp_culture_48h]),
    n_distinct(cohort_flags$patient_id[cohort_flags$any_resp_culture_48h & cohort_flags$any_positive_resp_culture_48h &
                                         cohort_flags$any_o2_ge2l_or_advanced_support_48h]),
    n_distinct(cohort_flags$patient_id[cohort_flags$final_cohort])
  )
)

definition_summary <- tibble(
  site_name = site_name,
  tables_path = tables_path,
  pathway = "adult ED -> ICU",
  time_window = "hospital admission through admission + 48h",
  respiratory_fluid_categories = paste(pulmonary_primary, collapse = ", "),
  positive_culture_definition = "method_category == culture and organism_group/category is not no_growth",
  oxygen_definition = "lpm_set >= 2 OR device_category in IMV/NIPPV/NIV/CPAP/HFNC/high flow NC/face mask/trach collar",
  antibiotic_definition = "administered med_category in antibacterial category list",
  n_final_hospitalizations = sum(cohort_flags$final_cohort, na.rm = TRUE),
  n_final_patients = n_distinct(cohort_flags$patient_id[cohort_flags$final_cohort]),
  n_final_with_cms_sepsis_antimicrobial = sum(
    cohort_flags$final_cohort & cohort_flags$any_cms_sepsis_qualifying_antimicrobial_48h,
    na.rm = TRUE
  )
)

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
summary_path <- file.path(out_dir, glue("ed_icu_early_resp_culture_cohort_count_{site_name}_{stamp}.csv"))
definition_path <- file.path(out_dir, glue("ed_icu_early_resp_culture_cohort_definition_{site_name}_{stamp}.csv"))

readr::write_csv(step_summary, summary_path)
readr::write_csv(definition_summary, definition_path)

message("Sequential cohort counts:")
print(step_summary, n = nrow(step_summary))
message("")
message("Definition summary:")
print(definition_summary)
message("Wrote count summary: ", summary_path)
message("Wrote definition summary: ", definition_path)
