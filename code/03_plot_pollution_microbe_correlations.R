# ================================================================================================
# CLIF Pollution-Microbiome | Correlation Tile Plot
# Purpose:
#   Plot organism_category correlations with PM2.5 and NO2 from the latest correlation output.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(glue)
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

out_dir <- file.path("output", "final")
fig_dir <- file.path("output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cor_path <- latest_file(file.path(out_dir, "pollution_microbe_correlations_*.csv"))
cor_dat <- readr::read_csv(cor_path, show_col_types = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

plot_dat <- cor_dat %>%
  filter(total_hospitalizations >= 25) %>%
  mutate(
    rank_score = pmax(
      abs(spearman_pm25_resp_cultured),
      abs(spearman_no2_resp_cultured),
      abs(spearman_pm25_positive),
      abs(spearman_no2_positive),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(rank_score), desc(total_hospitalizations)) %>%
  slice_head(n = 30) %>%
  select(
    organism_category,
    total_hospitalizations,
    `PM2.5 / respiratory-cultured` = spearman_pm25_resp_cultured,
    `NO2 / respiratory-cultured` = spearman_no2_resp_cultured,
    `PM2.5 / positive cultures` = spearman_pm25_positive,
    `NO2 / positive cultures` = spearman_no2_positive
  ) %>%
  pivot_longer(
    cols = starts_with(c("PM2.5", "NO2")),
    names_to = "contrast",
    values_to = "spearman_rho"
  ) %>%
  mutate(
    organism_label = glue("{pretty_organism(organism_category)} (n={total_hospitalizations})"),
    contrast = factor(
      contrast,
      levels = c(
        "PM2.5 / respiratory-cultured",
        "NO2 / respiratory-cultured",
        "PM2.5 / positive cultures",
        "NO2 / positive cultures"
      )
    )
  )

organism_order <- plot_dat %>%
  group_by(organism_label) %>%
  summarise(max_abs = max(abs(spearman_rho), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs) %>%
  pull(organism_label)

plot_dat <- plot_dat %>%
  mutate(organism_label = factor(organism_label, levels = organism_order))

p <- ggplot(plot_dat, aes(x = contrast, y = organism_label, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#B6423C",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,
    name = "Spearman rho"
  ) +
  labs(
    title = "County-Year Air Pollution vs Respiratory Culture Organisms",
    subtitle = "UCMC exploratory aggregate analysis; top 30 organisms by maximum absolute correlation",
    x = NULL,
    y = NULL,
    caption = glue("Correlation source: {basename(cor_path)}")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

png_path <- file.path(fig_dir, glue("pollution_microbe_correlation_tile_{stamp}.png"))
pdf_path <- file.path(fig_dir, glue("pollution_microbe_correlation_tile_{stamp}.pdf"))

ggsave(png_path, p, width = 10.5, height = 9, dpi = 300)
ggsave(pdf_path, p, width = 10.5, height = 9)

message("Wrote figure: ", png_path)
message("Wrote figure: ", pdf_path)
