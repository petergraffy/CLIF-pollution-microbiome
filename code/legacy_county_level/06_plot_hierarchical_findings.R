# ================================================================================================
# CLIF Pollution-Microbiome | Hierarchical Model Figures
# Purpose:
#   Create polished figures from patient-level hierarchical phenotype model outputs.
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

latest_file <- function(pattern) {
  hits <- Sys.glob(pattern)
  if (length(hits) == 0) stop("No files found for pattern: ", pattern)
  hits[which.max(file.info(hits)$mtime)]
}

pretty_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("pm25 mean", "PM2.5") %>%
    str_replace_all("no2 mean", "NO2") %>%
    str_replace_all("severe arf proxy", "Severe respiratory support") %>%
    str_to_sentence()
}

out_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

baseline_hits <- Sys.glob("output/final/hierarchical_pollution_microbe_phenotype_models_UCMC_all_counties_*.csv")
model_path <- if (length(baseline_hits) > 0) {
  baseline_hits[which.max(file.info(baseline_hits)$mtime)]
} else {
  legacy_hits <- Sys.glob("output/final/hierarchical_pollution_microbe_phenotype_models_UCMC_[0-9]*.csv")
  if (length(legacy_hits) == 0) latest_file("output/final/hierarchical_pollution_microbe_phenotype_models_UCMC_*.csv") else legacy_hits[which.max(file.info(legacy_hits)$mtime)]
}
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

model_results <- readr::read_csv(model_path, show_col_types = FALSE) %>%
  mutate(
    exposure_label = recode(exposure, pm25_mean = "PM2.5", no2_mean = "NO2"),
    phenotype_label = recode(
      phenotype,
      pneumonia_dx = "Pneumonia",
      sepsis_dx = "Sepsis",
      severe_arf_proxy = "Severe respiratory support"
    ),
    organism_label = pretty_label(organism_category),
    level_label = recode(phenotype_level, absent = "Phenotype absent", present = "Phenotype present"),
    signal = case_when(
      phenotype_level == "present" & odds_ratio_per_iqr > 1 & fdr_p_value < 0.10 ~ "FDR < 0.10",
      phenotype_level == "present" & odds_ratio_per_iqr > 1 & p_value < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "Other"
    )
  )

signal_orgs <- model_results %>%
  filter(phenotype_level == "present") %>%
  group_by(organism_category, organism_label) %>%
  summarise(best_p = min(p_value, na.rm = TRUE), max_or = max(odds_ratio_per_iqr, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_p, desc(max_or)) %>%
  slice_head(n = 12)

forest_dat <- model_results %>%
  filter(phenotype_level == "present", organism_category %in% signal_orgs$organism_category) %>%
  mutate(
    organism_label = factor(organism_label, levels = rev(signal_orgs$organism_label)),
    facet_label = paste(exposure_label, phenotype_label, sep = " / ")
  )

forest_dat <- forest_dat %>%
  mutate(
    odds_ratio_plot = case_when(
      exposure_label == "NO2" ~ pmin(pmax(odds_ratio_per_iqr, 0.80), 1.35),
      TRUE ~ pmin(pmax(odds_ratio_per_iqr, 0.45), 4.50)
    ),
    ci_low_plot = case_when(
      exposure_label == "NO2" ~ pmax(ci_low, 0.80),
      TRUE ~ pmax(ci_low, 0.45)
    ),
    ci_high_plot = case_when(
      exposure_label == "NO2" ~ pmin(ci_high, 1.35),
      TRUE ~ pmin(ci_high, 4.50)
    )
  )

forest <- ggplot(
  forest_dat,
  aes(x = odds_ratio_per_iqr, y = organism_label, color = signal)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(x = odds_ratio_plot, xmin = ci_low_plot, xmax = ci_high_plot), orientation = "y", width = 0.18, linewidth = 0.55) +
  geom_point(aes(x = odds_ratio_plot), size = 2.2) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4)) +
  scale_color_manual(
    values = c("FDR < 0.10" = "#B6423C", "Nominal p < 0.05" = "#D99A3D", "Other" = "grey45"),
    breaks = c("FDR < 0.10", "Nominal p < 0.05", "Other"),
    name = NULL
  ) +
  facet_grid(phenotype_label ~ exposure_label, scales = "free_x") +
  labs(
    title = "Pollution-Associated Respiratory Organism Detection by Clinical Phenotype",
    subtitle = "Mixed-effects logistic models with county random intercepts; odds ratios per IQR exposure increase",
    x = "Odds ratio per IQR increase in exposure",
    y = NULL,
    caption = glue("Model source: {basename(model_path)}")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.55),
    panel.spacing.x = unit(1.35, "lines"),
    strip.background = element_rect(fill = "grey92", color = "grey35", linewidth = 0.55),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

forest_png <- file.path(out_dir, glue("hierarchical_pollution_microbe_forest_{stamp}.png"))
forest_pdf <- file.path(out_dir, glue("hierarchical_pollution_microbe_forest_{stamp}.pdf"))
ggsave(forest_png, forest, width = 12.5, height = 9, dpi = 300)
ggsave(forest_pdf, forest, width = 12.5, height = 9)

heatmap_dat <- model_results %>%
  filter(phenotype_level == "present", organism_category %in% signal_orgs$organism_category) %>%
  mutate(
    organism_label = factor(organism_label, levels = rev(signal_orgs$organism_label)),
    contrast = factor(
      paste(exposure_label, phenotype_label, sep = " / "),
      levels = c(
        "PM2.5 / Pneumonia", "NO2 / Pneumonia",
        "PM2.5 / Sepsis", "NO2 / Sepsis",
        "PM2.5 / Severe respiratory support", "NO2 / Severe respiratory support"
      )
    ),
    tile_label = sprintf("%.2f", odds_ratio_per_iqr)
  )

heatmap <- ggplot(heatmap_dat, aes(x = contrast, y = organism_label, fill = log_or)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = tile_label), size = 3) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#B6423C",
    midpoint = 0,
    limits = c(log(0.5), log(2.5)),
    oob = scales::squish,
    name = "log(OR)"
  ) +
  labs(
    title = "Phenotype-Stratified Pollution Associations",
    subtitle = "Cells show odds ratio per IQR exposure increase",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

heatmap_png <- file.path(out_dir, glue("hierarchical_pollution_microbe_heatmap_{stamp}.png"))
heatmap_pdf <- file.path(out_dir, glue("hierarchical_pollution_microbe_heatmap_{stamp}.pdf"))
ggsave(heatmap_png, heatmap, width = 11, height = 7.5, dpi = 300)
ggsave(heatmap_pdf, heatmap, width = 11, height = 7.5)

interaction_dat <- model_results %>%
  filter(phenotype_level == "present", !is.na(interaction_or)) %>%
  group_by(organism_category, organism_label) %>%
  summarise(best_interaction_p = min(interaction_p_value, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_interaction_p) %>%
  slice_head(n = 12) %>%
  select(organism_category, organism_label) %>%
  inner_join(
    model_results %>%
      filter(phenotype_level == "present", !is.na(interaction_or)),
    by = c("organism_category", "organism_label")
  ) %>%
  mutate(
    organism_label = fct_reorder(organism_label, -interaction_p_value, .fun = min, na.rm = TRUE),
    contrast = factor(
      paste(exposure_label, phenotype_label, sep = " / "),
      levels = c(
        "PM2.5 / Pneumonia", "NO2 / Pneumonia",
        "PM2.5 / Sepsis", "NO2 / Sepsis",
        "PM2.5 / Severe respiratory support", "NO2 / Severe respiratory support"
      )
    )
  )

interaction_plot <- ggplot(interaction_dat, aes(x = interaction_or, y = organism_label, color = interaction_p_value < 0.05)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_point(size = 2) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  scale_color_manual(values = c("TRUE" = "#B6423C", "FALSE" = "grey45"), labels = c("FALSE" = "p >= 0.05", "TRUE" = "p < 0.05"), name = "Interaction") +
  facet_wrap(~ contrast, ncol = 2) +
  labs(
    title = "Pollution-by-Phenotype Interaction Signals",
    subtitle = "Interaction OR > 1 means the exposure-organism association is stronger when the phenotype is present",
    x = "Interaction odds ratio",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

interaction_png <- file.path(out_dir, glue("hierarchical_pollution_microbe_interactions_{stamp}.png"))
interaction_pdf <- file.path(out_dir, glue("hierarchical_pollution_microbe_interactions_{stamp}.pdf"))
ggsave(interaction_png, interaction_plot, width = 12, height = 9, dpi = 300)
ggsave(interaction_pdf, interaction_plot, width = 12, height = 9)

message("Wrote figures:")
message("  ", forest_png)
message("  ", heatmap_png)
message("  ", interaction_png)
