# ================================================================================================
# PheWAS-Style Plots for Positive Lung Culture Organism Groups vs Prior-Year Pollution
#
# Cohort:
#   Hospitalizations with any positive pulmonary culture in the CLIF database.
#
# Models:
#   For each organism_group, fit present vs absent among positive pulmonary culture
#   hospitalizations. Exposures are prior-year ZCTA PM2.5 and NO2, scaled as:
#     - PM2.5: 5 ug/m3
#     - NO2: 10 ppb
#
# Plot:
#   One lollipop-style PheWAS plot per exposure.
#   Y-axis is signed -log10(p): positive if OR > 1, negative if OR < 1.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(glue)
  library(readr)
  library(stringr)
  library(tidyr)
})

latest_file <- function(pattern, dir = file.path("output", "final")) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matched pattern: ", pattern)
  files[which.max(file.info(files)$mtime)]
}

latest_culture_exposure_file <- function(dir = file.path("output", "final")) {
  files <- list.files(
    dir,
    pattern = "^positive_lung_cultures_prior_year_pollution_.*\\.csv$",
    full.names = TRUE
  )
  files <- files[!str_detect(basename(files), "_coverage_|_organism_summary_|_year_summary_")]
  if (length(files) == 0) stop("No positive lung culture exposure exports found.")
  files[which.max(file.info(files)$mtime)]
}

clean_covariate <- function(x, missing = "unknown") {
  x <- str_to_lower(str_trim(as.character(x)))
  x <- if_else(is.na(x) | x == "", missing, x)
  factor(x)
}

fit_group_model <- function(data, group_name, exposure) {
  exposure_unit <- case_when(
    exposure == "pm25_prior_year" ~ 5,
    exposure == "no2_prior_year" ~ 10,
    TRUE ~ NA_real_
  )
  if (!is.finite(exposure_unit) || exposure_unit <= 0) stop("No exposure unit configured for: ", exposure)

  dat <- data %>%
    mutate(
      outcome = as.integer(organism_group == group_name),
      exposure_value = .data[[exposure]]
    ) %>%
    filter(!is.na(exposure_value), !is.na(outcome), !is.na(admission_year))

  n_events <- sum(dat$outcome == 1, na.rm = TRUE)
  n_nonevents <- sum(dat$outcome == 0, na.rm = TRUE)
  if (n_events < MIN_GROUP_DETECTIONS || n_nonevents < MIN_GROUP_DETECTIONS) return(NULL)

  dat <- dat %>%
    mutate(
      exposure_scaled = exposure_value / exposure_unit,
      admission_year = factor(admission_year),
      sex_category = clean_covariate(sex_category),
      race_category = fct_lump_min(clean_covariate(race_category), min = 50),
      ethnicity_category = clean_covariate(ethnicity_category)
    )

  fit <- tryCatch(
    suppressWarnings(glm(
      outcome ~ exposure_scaled + age_at_admission + sex_category +
        race_category + ethnicity_category + admission_year,
      family = binomial(),
      data = dat
    )),
    error = function(e) {
      warning("Model failed for ", group_name, " / ", exposure, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  coef_tab <- summary(fit)$coefficients
  if (!"exposure_scaled" %in% rownames(coef_tab)) return(NULL)

  beta <- coef_tab["exposure_scaled", "Estimate"]
  se <- coef_tab["exposure_scaled", "Std. Error"]
  p_value <- coef_tab["exposure_scaled", "Pr(>|z|)"]

  tibble(
    organism_group = group_name,
    exposure = exposure,
    n_hospitalizations = nrow(dat),
    n_events = n_events,
    n_controls = n_nonevents,
    exposure_unit = exposure_unit,
    exposure_unit_label = if_else(exposure == "pm25_prior_year", "5 ug/m3", "10 ppb"),
    odds_ratio_per_unit = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    p_value = p_value
  )
}

make_phewas_plot <- function(model_results, exposure_name, exposure_label, out_path_base) {
  plot_dat <- model_results %>%
    filter(exposure == exposure_name, !is.na(p_value), p_value > 0) %>%
    mutate(
      signed_log10_p = sign(log(odds_ratio_per_unit)) * -log10(p_value),
      direction = if_else(odds_ratio_per_unit >= 1, "Higher odds", "Lower odds"),
      organism_group_label = str_replace_all(organism_group, "_", " "),
      organism_group_label = fct_reorder(organism_group_label, signed_log10_p)
    )
  p <- ggplot(plot_dat, aes(x = organism_group_label, y = signed_log10_p, color = direction)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey45") +
    geom_hline(yintercept = c(-log10(0.05), log10(0.05)), linewidth = 0.3, linetype = "dashed", color = "grey65") +
    geom_segment(aes(xend = organism_group_label, y = 0, yend = signed_log10_p), linewidth = 0.45, alpha = 0.75) +
    geom_point(aes(size = n_events), alpha = 0.9) +
    scale_color_manual(values = c("Higher odds" = "#b23a48", "Lower odds" = "#2f6f8f")) +
    scale_size_continuous(range = c(1.8, 6.5), breaks = c(25, 100, 500, 1000, 2000), name = "Hospitalizations") +
    coord_flip() +
    labs(
      x = "Organism group",
      y = "Signed -log10(p)",
      color = NULL,
      title = glue("Positive lung culture organism groups vs prior-year {exposure_label}"),
      subtitle = "Reference lines: solid gray = OR 1; dashed gray = nominal p 0.05"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "plain", size = 13),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.title.x = element_text(margin = margin(t = 8))
    )

  ggsave(glue("{out_path_base}.png"), p, width = 9, height = 8.5, dpi = 300)
  ggsave(glue("{out_path_base}.jpg"), p, width = 9, height = 8.5, dpi = 300)
  ggsave(glue("{out_path_base}.pdf"), p, width = 9, height = 8.5)
  invisible(p)
}

MIN_GROUP_DETECTIONS <- as.integer(Sys.getenv("MIN_GROUP_DETECTIONS", unset = "10"))
input_path <- Sys.getenv(
  "POSITIVE_LUNG_CULTURE_EXPOSURE_PATH",
  unset = latest_culture_exposure_file()
)

message("Using culture exposure file: ", input_path)
message("Minimum group detections: ", MIN_GROUP_DETECTIONS)

culture_exposure <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    organism_group = str_to_lower(str_trim(as.character(organism_group))),
    admission_year = as.integer(admission_year),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission)),
    pm25_prior_year = suppressWarnings(as.numeric(pm25_prior_year)),
    no2_prior_year = suppressWarnings(as.numeric(no2_prior_year))
  ) %>%
  filter(!is.na(hospitalization_id), !is.na(organism_group), organism_group != "")

analysis_dat <- culture_exposure %>%
  distinct(
    hospitalization_id,
    organism_group,
    .keep_all = TRUE
  ) %>%
  select(
    hospitalization_id,
    patient_id,
    organism_group,
    admission_year,
    age_at_admission,
    sex_category,
    race_category,
    ethnicity_category,
    pm25_prior_year,
    no2_prior_year
  )

top_groups <- analysis_dat %>%
  distinct(hospitalization_id, organism_group) %>%
  count(organism_group, sort = TRUE, name = "n_events") %>%
  filter(n_events >= MIN_GROUP_DETECTIONS) %>%
  pull(organism_group)

message("Modeling ", length(top_groups), " organism groups")

modeled_group_counts <- analysis_dat %>%
  distinct(hospitalization_id, patient_id, organism_group) %>%
  filter(organism_group %in% top_groups) %>%
  group_by(organism_group) %>%
  summarise(
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_hospitalizations), organism_group)

model_grid <- tidyr::expand_grid(
  organism_group = top_groups,
  exposure = c("pm25_prior_year", "no2_prior_year")
)

model_results <- purrr::pmap_dfr(
  model_grid,
  ~ fit_group_model(analysis_dat, group_name = ..1, exposure = ..2)
) %>%
  group_by(exposure) %>%
  mutate(
    fdr_p_value = p.adjust(p_value, method = "BH"),
    signed_log10_p = sign(log(odds_ratio_per_unit)) * -log10(p_value)
  ) %>%
  ungroup() %>%
  arrange(exposure, p_value)

out_dir <- file.path("output", "final")
fig_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_path <- file.path(out_dir, glue("positive_lung_culture_group_prior_year_pollution_models_{stamp}.csv"))
summary_path <- file.path(out_dir, glue("positive_lung_culture_group_phewas_site_summary_{stamp}.csv"))
modeled_counts_path <- file.path(out_dir, glue("positive_lung_culture_group_phewas_modeled_groups_{stamp}.csv"))
readr::write_csv(model_results, model_path)

site_summary <- tibble(
  input_path = input_path,
  minimum_group_detections = MIN_GROUP_DETECTIONS,
  n_positive_lung_culture_rows = nrow(culture_exposure),
  n_positive_lung_culture_hospitalizations = n_distinct(culture_exposure$hospitalization_id),
  n_positive_lung_culture_patients = n_distinct(culture_exposure$patient_id),
  n_total_organism_groups = n_distinct(culture_exposure$organism_group),
  n_modeled_organism_groups = length(top_groups),
  n_model_rows = nrow(model_results),
  n_pm25_models = sum(model_results$exposure == "pm25_prior_year", na.rm = TRUE),
  n_no2_models = sum(model_results$exposure == "no2_prior_year", na.rm = TRUE)
)

readr::write_csv(site_summary, summary_path)
readr::write_csv(modeled_group_counts, modeled_counts_path)

pm25_plot_base <- file.path(fig_dir, glue("positive_lung_culture_group_phewas_pm25_prior_year_{stamp}"))
no2_plot_base <- file.path(fig_dir, glue("positive_lung_culture_group_phewas_no2_prior_year_{stamp}"))

make_phewas_plot(model_results, "pm25_prior_year", "PM2.5", pm25_plot_base)
make_phewas_plot(model_results, "no2_prior_year", "NO2", no2_plot_base)

message("Top associations by nominal p-value:")
print(
  model_results %>%
    select(organism_group, exposure, n_events, exposure_unit_label, odds_ratio_per_unit, ci_low, ci_high, p_value, fdr_p_value, signed_log10_p) %>%
    arrange(p_value) %>%
    head(20),
  n = 20
)
message("Wrote model results: ", model_path)
message("Wrote site summary: ", summary_path)
message("Wrote modeled group counts: ", modeled_counts_path)
message("Wrote PM2.5 plot: ", pm25_plot_base, ".png")
message("Wrote NO2 plot: ", no2_plot_base, ".png")
message("Wrote PM2.5 JPG: ", pm25_plot_base, ".jpg")
message("Wrote NO2 JPG: ", no2_plot_base, ".jpg")
