# ================================================================================================
# PBF grant figure | Positive lung culture organism-wide association screen
# Purpose:
#   Create a grant-ready prior-year NO2 organism PheWAS plot with larger text,
#   quieter non-significant points, and no caption/subtitle.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(glue)
  library(readr)
  library(stringr)
})

model_path <- Sys.getenv(
  "PBF_ORGANISM_PHEWAS_MODEL_PATH",
  unset = "output/final/positive_lung_culture_organism_prior_year_pollution_models_20260727_150720.csv"
)

if (!file.exists(model_path)) {
  stop("Model file not found: ", model_path)
}

pretty_taxon <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("baumanii", "baumannii") %>%
    str_replace_all("intercelluare", "intracellulare")
}

out_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_dat <- readr::read_csv(model_path, show_col_types = FALSE) %>%
  filter(exposure == "no2_prior_year", !is.na(p_value), p_value > 0) %>%
  mutate(
    neg_log10_p = -log10(p_value),
    organism_label = pretty_taxon(organism_category),
    group_label = pretty_taxon(organism_group),
    direction = if_else(odds_ratio_per_unit >= 1, "Higher odds", "Lower odds"),
    signal = case_when(
      fdr_p_value < 0.10 ~ "FDR < 0.10",
      p_value < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(organism_group, desc(neg_log10_p), organism_category) %>%
  group_by(organism_group) %>%
  mutate(group_rank = row_number()) %>%
  ungroup() %>%
  arrange(organism_group, group_rank) %>%
  mutate(x_pos = row_number())

group_bounds <- plot_dat %>%
  group_by(organism_group, group_label) %>%
  summarise(
    xmin = min(x_pos) - 0.5,
    xmax = max(x_pos) + 0.5,
    xcenter = mean(x_pos),
    .groups = "drop"
  ) %>%
  mutate(
    group_index = row_number(),
    shade = group_index %% 2 == 1,
    group_label_wrapped = str_wrap(group_label, width = 13)
  )

threshold_p05 <- -log10(0.05)
threshold_fdr <- plot_dat %>%
  filter(fdr_p_value < 0.10) %>%
  summarise(y = min(neg_log10_p, na.rm = TRUE)) %>%
  pull(y)
threshold_fdr <- ifelse(length(threshold_fdr) == 0 || !is.finite(threshold_fdr), NA_real_, threshold_fdr)

y_max <- max(plot_dat$neg_log10_p, threshold_p05, threshold_fdr, na.rm = TRUE)
x_label <- max(plot_dat$x_pos, na.rm = TRUE)

label_dat <- plot_dat %>%
  filter(p_value < 0.01) %>%
  mutate(
    label_nudge_y = if_else(neg_log10_p > 4.2, 0.10, 0.05)
  )

phewas <- ggplot(plot_dat, aes(x = x_pos, y = neg_log10_p)) +
  geom_rect(
    data = group_bounds %>% filter(shade),
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey94",
    alpha = 0.7
  ) +
  geom_hline(yintercept = threshold_p05, linewidth = 0.65, color = "#D95F5F", alpha = 0.8) +
  {if (!is.na(threshold_fdr)) geom_hline(
    yintercept = threshold_fdr,
    linewidth = 0.65,
    linetype = "dashed",
    color = "#8B1E3F",
    alpha = 0.85
  )} +
  annotate(
    "text",
    x = x_label,
    y = threshold_p05,
    label = "Nominal p = 0.05",
    hjust = 1,
    vjust = -0.45,
    size = 5.2,
    color = "#D95F5F"
  ) +
  {if (!is.na(threshold_fdr)) annotate(
    "text",
    x = x_label,
    y = threshold_fdr,
    label = "FDR < 0.10 threshold",
    hjust = 1,
    vjust = -0.45,
    size = 5.2,
    color = "#8B1E3F"
  )} +
  geom_point(
    aes(color = signal, shape = direction, size = n_events),
    stroke = 1.05,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = label_dat,
    aes(label = organism_label),
    size = 4.5,
    color = "black",
    seed = 1729,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.size = 0.35,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    force = 3,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = group_bounds$xcenter,
    labels = NULL,
    expand = expansion(mult = c(0.012, 0.025))
  ) +
  scale_y_continuous(
    limits = c(0, y_max * 1.12),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_color_manual(
    values = c(
      "FDR < 0.10" = "#B6423C",
      "Nominal p < 0.05" = "#D99A3D",
      "Other" = "grey60"
    ),
    breaks = c("FDR < 0.10", "Nominal p < 0.05", "Other"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Higher odds" = 24, "Lower odds" = 25),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(2.0, 7.0),
    breaks = c(25, 100, 500, 1000),
    name = "Detections"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4.5)),
    shape = guide_legend(order = 2, override.aes = list(size = 4.5)),
    size = guide_legend(order = 3)
  ) +
  labs(
    title = "Prior-Year NO2 Exposure and Respiratory Organism Detection",
    x = "Respiratory organisms grouped by taxonomy",
    y = expression(-log[10](p))
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.55),
    plot.title = element_text(face = "bold", size = 22, color = "black", margin = margin(b = 8)),
    axis.title.x = element_text(size = 17, color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 17, color = "black", margin = margin(r = 8)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 13, color = "grey25"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.key.width = unit(1.0, "lines"),
    plot.margin = margin(t = 12, r = 16, b = 8, l = 14)
  )

png_path <- file.path(out_dir, "pbf_positive_lung_organism_phewas_no2_grant.png")
pdf_path <- file.path(out_dir, "pbf_positive_lung_organism_phewas_no2_grant.pdf")

ggsave(png_path, phewas, width = 15.5, height = 9.0, dpi = 450, bg = "white")
ggsave(pdf_path, phewas, width = 15.5, height = 9.0, bg = "white")

message("Wrote grant-ready organism PheWAS figure:")
message("  ", normalizePath(png_path))
message("  ", normalizePath(pdf_path))
