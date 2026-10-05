#!/usr/bin/env Rscript

# Grant-ready ZCTA maps of study-period air pollution exposure surfaces.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(patchwork)
  library(readr)
  library(sf)
  library(tibble)
  library(viridis)
})

sf_use_s2(FALSE)

start_year <- as.integer(Sys.getenv("POLLUTION_MAP_START_YEAR", "2018"))
end_year <- as.integer(Sys.getenv("POLLUTION_MAP_END_YEAR", "2025"))

out_dir <- file.path("output", "figures", "study_period_pollution_maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zcta_path <- Sys.getenv(
  "ZCTA_SHAPEFILE",
  unset = "/Users/saborpete/Desktop/Peter/Postdoc/environment_transplant_survival/data/cache/cb_2020_us_zcta520_500k/cb_2020_us_zcta520_500k.shp"
)

pollutants <- tribble(
  ~pollutant, ~file, ~value_col, ~label, ~legend, ~palette, ~direction,
  "pm25",
  "data/exposome_zcta/air_pollution_zcta_pm25_monthly_2005_2023.parquet",
  "pm25_ug_m3",
  "Fine particulate matter",
  expression("Fine particulate matter ("*mu*"g/m"^3*")"),
  "magma",
  1,
  "no2",
  "data/exposome_zcta/air_pollution_zcta_no2_annual_2005_2025.parquet",
  "no2",
  "Nitrogen dioxide",
  "Nitrogen dioxide (ppb)",
  "inferno",
  1,
  "o3",
  "data/exposome_zcta/air_pollution_zcta_o3_monthly_2005_2023.parquet",
  "o3_ppb",
  "Ozone",
  "Ozone (ppb)",
  "viridis",
  1
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

read_pollutant_mean <- function(row) {
  path <- row$file
  value_col <- row$value_col
  if (!file.exists(path)) stop("Missing pollutant file: ", path, call. = FALSE)

  ds <- arrow::open_dataset(path)
  available <- ds %>%
    summarise(first_year = min(year), last_year = max(year)) %>%
    collect()

  use_start <- max(start_year, available$first_year)
  use_end <- min(end_year, available$last_year)
  if (use_start > use_end) {
    stop("No years available for ", row$pollutant, " in requested window.", call. = FALSE)
  }

  log_msg("Summarising ", row$label, " for ", use_start, "-", use_end)
  dat <- ds %>%
    filter(year >= use_start, year <= use_end) %>%
    select(zip, year, value = all_of(value_col)) %>%
    group_by(zip) %>%
    summarise(
      value = mean(value, na.rm = TRUE),
      n_records = n(),
      first_year = min(year),
      last_year = max(year),
      .groups = "drop"
    ) %>%
    collect() %>%
    mutate(
      pollutant = row$pollutant,
      label = row$label,
      units = as.character(row$legend),
      requested_start_year = start_year,
      requested_end_year = end_year,
      available_start_year = available$first_year,
      available_end_year = available$last_year,
      mapped_start_year = use_start,
      mapped_end_year = use_end
    )

  dat
}

if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

log_msg("Reading ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
conus_bbox <- st_transform(conus_bbox, st_crs(zcta))
zcta <- suppressWarnings(st_crop(zcta, conus_bbox))
zcta <- st_transform(zcta, 5070)

pollution_means <- bind_rows(lapply(seq_len(nrow(pollutants)), function(i) read_pollutant_mean(pollutants[i, ])))

write_csv(
  pollution_means,
  file.path(out_dir, glue("zcta_pollutant_means_{start_year}_{end_year}_available_years.csv.gz"))
)

summary_tbl <- pollution_means %>%
  group_by(pollutant, label, mapped_start_year, mapped_end_year, available_start_year, available_end_year) %>%
  summarise(
    zctas = n_distinct(zip),
    min = min(value, na.rm = TRUE),
    p01 = quantile(value, 0.01, na.rm = TRUE),
    p05 = quantile(value, 0.05, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    p95 = quantile(value, 0.95, na.rm = TRUE),
    p99 = quantile(value, 0.99, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_tbl, file.path(out_dir, glue("study_period_pollution_map_summary_{start_year}_{end_year}.csv")))

plot_one <- function(pollutant_id) {
  meta <- pollutants %>% filter(pollutant == pollutant_id) %>% slice(1)
  dat <- pollution_means %>% filter(pollutant == pollutant_id)
  yrs <- dat %>% summarise(first = min(mapped_start_year), last = max(mapped_end_year))
  lims <- quantile(dat$value, c(0.02, 0.98), na.rm = TRUE)
  map_data <- merge(zcta, dat %>% select(zip, value), by = "zip", all.x = TRUE)
  map_data <- map_data %>% mutate(value_plot = pmin(pmax(value, lims[[1]]), lims[[2]]))

  missing_polygons <- sum(is.na(map_data$value))
  write_csv(
    tibble(
      pollutant = pollutant_id,
      map_polygons = nrow(map_data),
      missing_polygons = missing_polygons,
      mapped_start_year = yrs$first,
      mapped_end_year = yrs$last,
      color_scale_p02 = unname(lims[[1]]),
      color_scale_p98 = unname(lims[[2]])
    ),
    file.path(out_dir, glue("{pollutant_id}_map_manifest_{start_year}_{end_year}.csv"))
  )

  title_year <- if (yrs$first == yrs$last) as.character(yrs$first) else glue("{yrs$first}-{yrs$last}")

  ggplot(map_data) +
    geom_sf(aes(fill = value_plot), color = NA) +
    scale_fill_viridis(
      option = meta$palette,
      direction = meta$direction,
      name = meta$legend[[1]],
      limits = lims,
      na.value = "grey92",
      guide = guide_colorbar(
        barheight = unit(42, "mm"),
        barwidth = unit(4.5, "mm"),
        ticks.colour = "grey30",
        frame.colour = "grey60"
      )
    ) +
    coord_sf(datum = NA) +
    labs(title = NULL, subtitle = NULL) +
    theme_void(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_text(size = 11, color = "black"),
      legend.text = element_text(size = 10, color = "black"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

log_msg("Rendering individual maps")
plots <- setNames(lapply(pollutants$pollutant, plot_one), pollutants$pollutant)

for (nm in names(plots)) {
  png_path <- file.path(out_dir, glue("{nm}_study_period_zcta_map_{start_year}_{end_year}.png"))
  pdf_path <- file.path(out_dir, glue("{nm}_study_period_zcta_map_{start_year}_{end_year}.pdf"))
  ggsave(png_path, plots[[nm]], width = 10.5, height = 6.6, dpi = 450, bg = "white")
  ggsave(pdf_path, plots[[nm]], width = 10.5, height = 6.6, bg = "white")
}

panel <- (plots$pm25 | plots$no2 | plots$o3) +
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    theme = theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  )

ggsave(
  file.path(out_dir, glue("study_period_pollution_zcta_map_panel_{start_year}_{end_year}.png")),
  panel,
  width = 18,
  height = 4.8,
  dpi = 450,
  bg = "white"
)
ggsave(
  file.path(out_dir, glue("study_period_pollution_zcta_map_panel_{start_year}_{end_year}.pdf")),
  panel,
  width = 18,
  height = 4.8,
  bg = "white"
)

manifest <- summary_tbl %>%
  transmute(
    pollutant,
    map_label = label,
    requested_study_period = glue("{start_year}-{end_year}"),
    mapped_years = glue("{mapped_start_year}-{mapped_end_year}"),
    available_years = glue("{available_start_year}-{available_end_year}"),
    zctas,
    median,
    p05,
    p95
  )
write_csv(manifest, file.path(out_dir, glue("study_period_pollution_maps_manifest_{start_year}_{end_year}.csv")))

log_msg("Wrote pollutant maps to ", normalizePath(out_dir, winslash = "/"))
