# ================================================================================================
# PBF grant figure | Pollution-associated respiratory organism detection
# Purpose:
#   Create a grant-ready forest plot with larger text and no subtitle/caption.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(glue)
  library(readr)
  library(stringr)
})

pretty_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("baumanii", "baumannii") %>%
    str_to_sentence()
}

model_path <- Sys.getenv(
  "PBF_MODEL_PATH",
  unset = "output/final/hierarchical_pollution_microbe_phenotype_models_UCMC_20260522_134554.csv"
)

if (!file.exists(model_path)) {
  stop("Model file not found: ", model_path)
}

out_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_results <- readr::read_csv(model_path, show_col_types = FALSE) %>%
  mutate(
    exposure_label = recode(exposure, pm25_mean = "PM2.5", no2_mean = "NO2"),
    phenotype_label = recode(
      phenotype,
      pneumonia_dx = "Pneumonia",
      sepsis_dx = "Sepsis",
      severe_arf_proxy = "Severe resp. support"
    ),
    organism_label = pretty_label(organism_category),
    signal = case_when(
      phenotype_level == "present" & odds_ratio_per_iqr > 1 & fdr_p_value < 0.10 ~ "FDR < 0.10",
      phenotype_level == "present" & odds_ratio_per_iqr > 1 & p_value < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "Other"
    )
  )

signal_orgs <- model_results %>%
  filter(phenotype_level == "present") %>%
  group_by(organism_category, organism_label) %>%
  summarise(
    best_p = min(p_value, na.rm = TRUE),
    max_or = max(odds_ratio_per_iqr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(best_p, desc(max_or)) %>%
  slice_head(n = 12)

forest_dat <- model_results %>%
  filter(phenotype_level == "present", organism_category %in% signal_orgs$organism_category) %>%
  mutate(
    organism_label = factor(organism_label, levels = rev(signal_orgs$organism_label)),
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
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.55) +
  geom_errorbar(
    aes(x = odds_ratio_plot, xmin = ci_low_plot, xmax = ci_high_plot),
    orientation = "y",
    width = 0.18,
    linewidth = 0.75
  ) +
  geom_point(aes(x = odds_ratio_plot), size = 3.0) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4),
    labels = c("0.50", "0.75", "1.00", "1.25", "1.50", "2.00", "3.00", "4.00")
  ) +
  scale_color_manual(
    values = c(
      "FDR < 0.10" = "#B6423C",
      "Nominal p < 0.05" = "#D99A3D",
      "Other" = "grey45"
    ),
    breaks = c("FDR < 0.10", "Nominal p < 0.05", "Other"),
    name = NULL
  ) +
  facet_grid(phenotype_label ~ exposure_label, scales = "free_x") +
  labs(
    title = "Pollution-Associated Respiratory Organism Detection by Clinical Phenotype",
    x = "Odds ratio per IQR increase in exposure",
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.65),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(0.65, "lines"),
    strip.background = element_rect(fill = "grey92", color = "grey35", linewidth = 0.65),
    strip.text = element_text(face = "bold", size = 13, color = "black"),
    strip.text.y = element_text(angle = 270, margin = margin(l = 6, r = 6)),
    plot.title = element_text(face = "bold", size = 18, color = "black", margin = margin(b = 9)),
    axis.title.x = element_text(size = 15, color = "black", margin = margin(t = 7)),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "grey20"),
    legend.position = "bottom",
    legend.text = element_text(size = 13, color = "black"),
    legend.key.width = unit(1.2, "lines"),
    legend.margin = margin(t = 2),
    plot.margin = margin(t = 10, r = 12, b = 8, l = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.4, linewidth = 0.9)))

png_path <- file.path(out_dir, "pbf_pollution_microbe_forest_grant.png")
pdf_path <- file.path(out_dir, "pbf_pollution_microbe_forest_grant.pdf")

ggsave(png_path, forest, width = 12.8, height = 8.4, dpi = 450, bg = "white")
ggsave(pdf_path, forest, width = 12.8, height = 8.4, bg = "white")

message("Wrote grant-ready figure:")
message("  ", normalizePath(png_path))
message("  ", normalizePath(pdf_path))
