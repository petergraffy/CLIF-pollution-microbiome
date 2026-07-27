# ================================================================================================
# PheWAS-Style Positive Lung Culture Organism Screen
#
# Cohort:
#   Hospitalizations with any positive pulmonary culture in the CLIF database.
#
# Models:
#   For each organism_category, fit organism present vs absent among positive pulmonary culture
#   hospitalizations. Exposures are prior-year ZCTA PM2.5 and NO2, scaled as:
#     - PM2.5: 5 ug/m3
#     - NO2: 10 ppb
#
# Plot:
#   PheWAS-style Manhattan plots with individual organisms arranged under organism_group.
#   Y-axis is -log10(p). Point shape shows direction of association.
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

fit_organism_model <- function(data, organism_name, exposure) {
  exposure_unit <- case_when(
    exposure == "pm25_prior_year" ~ 5,
    exposure == "no2_prior_year" ~ 10,
    TRUE ~ NA_real_
  )
  if (!is.finite(exposure_unit) || exposure_unit <= 0) stop("No exposure unit configured for: ", exposure)

  dat <- data %>%
    mutate(
      outcome = as.integer(organism_category == organism_name),
      exposure_value = .data[[exposure]]
    ) %>%
    filter(!is.na(exposure_value), !is.na(outcome), !is.na(admission_year))

  n_events <- sum(dat$outcome == 1, na.rm = TRUE)
  n_nonevents <- sum(dat$outcome == 0, na.rm = TRUE)
  if (n_events < MIN_ORGANISM_DETECTIONS || n_nonevents < MIN_ORGANISM_DETECTIONS) return(NULL)

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
      warning("Model failed for ", organism_name, " / ", exposure, ": ", conditionMessage(e))
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
    organism_category = organism_name,
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
    arrange(organism_group, desc(neg_log10_p), organism_category) %>%
    group_by(organism_group) %>%
    mutate(group_rank = row_number()) %>%
    ungroup() %>%
    arrange(organism_group, group_rank) %>%
    mutate(
      x_pos = row_number(),
      organism_label = str_replace_all(organism_category, "_", " "),
      group_label = str_replace_all(organism_group, "_", " "),
      direction = if_else(odds_ratio_per_unit >= 1, "Higher odds", "Lower odds"),
      label_hit = fdr_p_value < 0.10 | p_value < 0.01
    )

  group_bounds <- plot_dat %>%
    group_by(organism_group, group_label) %>%
    summarise(
      xmin = min(x_pos) - 0.5,
      xmax = max(x_pos) + 0.5,
      xcenter = mean(x_pos),
      .groups = "drop"
    ) %>%
    mutate(group_index = row_number(), shade = group_index %% 2 == 1)

  threshold_p05 <- -log10(0.05)
  threshold_fdr <- plot_dat %>%
    filter(fdr_p_value < 0.10) %>%
    summarise(y = min(neg_log10_p, na.rm = TRUE)) %>%
    pull(y)
  threshold_fdr <- ifelse(length(threshold_fdr) == 0 || !is.finite(threshold_fdr), NA_real_, threshold_fdr)
  x_label <- max(plot_dat$x_pos, na.rm = TRUE)
  y_max <- max(plot_dat$neg_log10_p, threshold_p05, threshold_fdr, na.rm = TRUE)

  p <- ggplot(plot_dat, aes(x = x_pos, y = neg_log10_p)) +
    geom_rect(
      data = group_bounds %>% filter(shade),
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey93",
      alpha = 0.55
    ) +
    geom_hline(yintercept = threshold_p05, linewidth = 0.35, color = "#d95f5f", alpha = 0.75) +
    {if (!is.na(threshold_fdr)) geom_hline(yintercept = threshold_fdr, linewidth = 0.35, linetype = "dashed", color = "#8b1e3f", alpha = 0.75)} +
    annotate(
      "text",
      x = x_label,
      y = threshold_p05,
      label = "Nominal p = 0.05",
      hjust = 1,
      vjust = -0.45,
      size = 3,
      color = "#d95f5f"
    ) +
    {if (!is.na(threshold_fdr)) annotate(
      "text",
      x = x_label,
      y = threshold_fdr,
      label = "FDR < 0.10 threshold",
      hjust = 1,
      vjust = -0.45,
      size = 3,
      color = "#8b1e3f"
    )} +
    geom_point(aes(color = organism_group, shape = direction, size = n_events), alpha = 0.9) +
    geom_text(
      data = plot_dat %>% filter(label_hit),
      aes(label = organism_label),
      size = 2.6,
      angle = 35,
      hjust = -0.05,
      vjust = 0.4,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = group_bounds$xcenter,
      labels = group_bounds$group_label,
      expand = expansion(mult = c(0.015, 0.025))
    ) +
    scale_y_continuous(limits = c(0, y_max * 1.08), expand = expansion(mult = c(0.02, 0.04))) +
    scale_shape_manual(values = c("Higher odds" = 24, "Lower odds" = 25)) +
    scale_size_continuous(range = c(1.4, 5.2), breaks = c(25, 100, 500, 1000, 2000), name = "Detections") +
    guides(color = "none", shape = guide_legend(title = NULL), size = guide_legend(title = "Detections")) +
    labs(
      x = "Organism group",
      y = expression(-log[10](p)),
      title = glue("Organism-wide association screen: prior-year {exposure_label}")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "plain", size = 13),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.3),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.position = "bottom"
    )

  ggsave(glue("{out_path_base}.png"), p, width = 13.5, height = 8.2, dpi = 300)
  ggsave(glue("{out_path_base}.jpg"), p, width = 13.5, height = 8.2, dpi = 300)
  ggsave(glue("{out_path_base}.pdf"), p, width = 13.5, height = 8.2)
  invisible(p)
}

MIN_ORGANISM_DETECTIONS <- as.integer(Sys.getenv("MIN_ORGANISM_DETECTIONS", unset = "10"))
input_path <- Sys.getenv(
  "POSITIVE_LUNG_CULTURE_EXPOSURE_PATH",
  unset = latest_culture_exposure_file()
)

message("Using culture exposure file: ", input_path)
message("Minimum organism detections: ", MIN_ORGANISM_DETECTIONS)

culture_exposure <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    organism_category = str_to_lower(str_trim(as.character(organism_category))),
    organism_group = str_to_lower(str_trim(as.character(organism_group))),
    admission_year = as.integer(admission_year),
    age_at_admission = suppressWarnings(as.numeric(age_at_admission)),
    pm25_prior_year = suppressWarnings(as.numeric(pm25_prior_year)),
    no2_prior_year = suppressWarnings(as.numeric(no2_prior_year))
  ) %>%
  filter(!is.na(hospitalization_id), !is.na(organism_category), organism_category != "")

organism_group_map <- culture_exposure %>%
  distinct(hospitalization_id, organism_category, organism_group) %>%
  count(organism_category, organism_group, sort = TRUE, name = "n") %>%
  group_by(organism_category) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(organism_category, organism_group)

analysis_dat <- culture_exposure %>%
  distinct(
    hospitalization_id,
    organism_category,
    .keep_all = TRUE
  ) %>%
  select(
    hospitalization_id,
    patient_id,
    organism_category,
    admission_year,
    age_at_admission,
    sex_category,
    race_category,
    ethnicity_category,
    pm25_prior_year,
    no2_prior_year
  )

top_organisms <- analysis_dat %>%
  distinct(hospitalization_id, organism_category) %>%
  count(organism_category, sort = TRUE, name = "n_events") %>%
  filter(n_events >= MIN_ORGANISM_DETECTIONS) %>%
  pull(organism_category)

message("Modeling ", length(top_organisms), " organisms")

modeled_organism_counts <- analysis_dat %>%
  distinct(hospitalization_id, patient_id, organism_category) %>%
  filter(organism_category %in% top_organisms) %>%
  left_join(organism_group_map, by = "organism_category") %>%
  mutate(organism_group = coalesce(organism_group, "ungrouped")) %>%
  group_by(organism_category, organism_group) %>%
  summarise(
    n_hospitalizations = n_distinct(hospitalization_id),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_hospitalizations), organism_group, organism_category)

model_grid <- tidyr::expand_grid(
  organism_category = top_organisms,
  exposure = c("pm25_prior_year", "no2_prior_year")
)

model_results <- purrr::pmap_dfr(
  model_grid,
  ~ fit_organism_model(analysis_dat, organism_name = ..1, exposure = ..2)
) %>%
  left_join(organism_group_map, by = "organism_category") %>%
  mutate(organism_group = coalesce(organism_group, "ungrouped")) %>%
  group_by(exposure) %>%
  mutate(
    fdr_p_value = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(p_value)
  ) %>%
  ungroup() %>%
  arrange(exposure, organism_group, p_value)

out_dir <- file.path("output", "final")
fig_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_path <- file.path(out_dir, glue("positive_lung_culture_organism_prior_year_pollution_models_{stamp}.csv"))
summary_path <- file.path(out_dir, glue("positive_lung_culture_organism_phewas_site_summary_{stamp}.csv"))
modeled_counts_path <- file.path(out_dir, glue("positive_lung_culture_organism_phewas_modeled_organisms_{stamp}.csv"))
readr::write_csv(model_results, model_path)

site_summary <- tibble(
  input_path = input_path,
  minimum_organism_detections = MIN_ORGANISM_DETECTIONS,
  n_positive_lung_culture_rows = nrow(culture_exposure),
  n_positive_lung_culture_hospitalizations = n_distinct(culture_exposure$hospitalization_id),
  n_positive_lung_culture_patients = n_distinct(culture_exposure$patient_id),
  n_total_organism_categories = n_distinct(culture_exposure$organism_category),
  n_total_organism_groups = n_distinct(culture_exposure$organism_group),
  n_modeled_organisms = length(top_organisms),
  n_modeled_organism_groups = n_distinct(modeled_organism_counts$organism_group),
  n_model_rows = nrow(model_results),
  n_pm25_models = sum(model_results$exposure == "pm25_prior_year", na.rm = TRUE),
  n_no2_models = sum(model_results$exposure == "no2_prior_year", na.rm = TRUE)
)

readr::write_csv(site_summary, summary_path)
readr::write_csv(modeled_organism_counts, modeled_counts_path)

pm25_plot_base <- file.path(fig_dir, glue("positive_lung_culture_organism_phewas_pm25_prior_year_{stamp}"))
no2_plot_base <- file.path(fig_dir, glue("positive_lung_culture_organism_phewas_no2_prior_year_{stamp}"))

make_phewas_plot(model_results, "pm25_prior_year", "PM2.5", pm25_plot_base)
make_phewas_plot(model_results, "no2_prior_year", "NO2", no2_plot_base)

message("Top associations by nominal p-value:")
print(
  model_results %>%
    select(organism_category, organism_group, exposure, n_events, exposure_unit_label, odds_ratio_per_unit, ci_low, ci_high, p_value, fdr_p_value) %>%
    arrange(p_value) %>%
    head(25),
  n = 25
)
message("Wrote model results: ", model_path)
message("Wrote site summary: ", summary_path)
message("Wrote modeled organism counts: ", modeled_counts_path)
message("Wrote PM2.5 plot: ", pm25_plot_base, ".png")
message("Wrote NO2 plot: ", no2_plot_base, ".png")
message("Wrote PM2.5 JPG: ", pm25_plot_base, ".jpg")
message("Wrote NO2 JPG: ", no2_plot_base, ".jpg")
