# ================================================================================================
# SHRF Pollution-Organism Forest Plot
# Purpose:
#   Plot organism-specific first-48h pulmonary culture associations with ZCTA pollution.
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

pretty_organism <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

model_path <- Sys.getenv(
  "SHRF_ORGANISM_MODEL_PATH",
  unset = latest_file("output/final/shrf_zcta_pollution_organism_models_*.csv")
)
out_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

model_results <- readr::read_csv(model_path, show_col_types = FALSE) %>%
  mutate(
    exposure_label = recode(
      exposure,
      pm25_annual = "PM2.5",
      o3_annual = "O3",
      no2_annual = "NO2",
      .default = exposure
    ),
    organism_label = pretty_organism(organism_category),
    event_label = glue("n={n_events}"),
    signal = case_when(
      odds_ratio_per_iqr > 1 & fdr_p_value < 0.10 ~ "FDR < 0.10",
      odds_ratio_per_iqr > 1 & p_value < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "Other"
    )
  )

organism_order <- model_results %>%
  group_by(organism_category, organism_label) %>%
  summarise(max_or = max(odds_ratio_per_iqr, na.rm = TRUE), best_p = min(p_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    max_or = if_else(is.finite(max_or), max_or, 0),
    best_p = if_else(is.finite(best_p), best_p, Inf)
  ) %>%
  arrange(desc(max_or), best_p) %>%
  pull(organism_label)

plot_dat <- tidyr::expand_grid(
  organism_category = unique(model_results$organism_category),
  exposure_label = factor(c("PM2.5", "O3", "NO2"), levels = c("PM2.5", "O3", "NO2"))
) %>%
  left_join(model_results %>% distinct(organism_category, organism_label), by = "organism_category") %>%
  left_join(
    model_results %>% select(-organism_label),
    by = c("organism_category", "exposure_label")
  ) %>%
  mutate(
    organism_label = factor(organism_label, levels = rev(unique(organism_order))),
    exposure_label = factor(exposure_label, levels = c("PM2.5", "O3", "NO2")),
    event_label = if_else(is.na(n_events), "not modeled", glue("n={n_events}")),
    signal = coalesce(signal, "Not modeled")
  )

legend_breaks <- c("FDR < 0.10", "Nominal p < 0.05", "Other")

modeled_dat <- plot_dat %>% filter(!is.na(odds_ratio_per_iqr))

event_x <- plot_dat %>%
  group_by(exposure_label) %>%
  summarise(
    label_x = max(ci_high, odds_ratio_per_iqr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label_x = if_else(is.finite(label_x), label_x, 1.5)
  )

plot_dat <- plot_dat %>% left_join(event_x, by = "exposure_label")

forest <- ggplot(plot_dat, aes(y = organism_label, color = signal)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_errorbar(data = modeled_dat, aes(x = odds_ratio_per_iqr, xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18, linewidth = 0.62) +
  geom_point(data = modeled_dat, aes(x = odds_ratio_per_iqr), size = 2.5) +
  geom_text(
    data = modeled_dat %>% left_join(event_x, by = "exposure_label"),
    aes(x = label_x, label = event_label),
    hjust = -0.08,
    color = "grey25",
    size = 2.9,
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3, 5),
    labels = c("0.5", "0.75", "1.0", "1.25", "1.5", "2.0", "3.0", "5.0")
  ) +
  scale_color_manual(
    values = c(
      "FDR < 0.10" = "#B6423C",
      "Nominal p < 0.05" = "#D9902F",
      "Other" = "#3F5968",
      "Not modeled" = "grey60"
    ),
    breaks = legend_breaks,
    name = NULL
  ) +
  facet_wrap(~ exposure_label, ncol = 3, scales = "free_x") +
  coord_cartesian(clip = "off") +
  labs(
    title = "Pollution and Early Pulmonary Organism Detection in Severe Hypoxemic Respiratory Failure",
    x = "Odds ratio per IQR exposure increase",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.55),
    panel.spacing.x = unit(2, "lines"),
    strip.background = element_rect(fill = "grey90", color = "grey35", linewidth = 0.55),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    plot.margin = margin(10, 42, 10, 10)
  )

png_path <- file.path(out_dir, glue("shrf_pollution_organism_forest_{stamp}.png"))
pdf_path <- file.path(out_dir, glue("shrf_pollution_organism_forest_{stamp}.pdf"))

ggsave(png_path, forest, width = 15, height = 7.5, dpi = 300)
ggsave(pdf_path, forest, width = 15, height = 7.5)

message("Wrote forest plot:")
message("  ", png_path)
message("  ", pdf_path)
