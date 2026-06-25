# ================================================================================================
# CLIF Pollution-Microbiome | Pollution-Organism Risk Models
# Purpose:
#   Fit aggregate binomial models testing whether higher county-year PM2.5 or NO2 corresponds to
#   higher odds of detecting specific respiratory culture organisms.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
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

fit_one_model <- function(df, organism, exposure, denominator_col, model_type) {
  dat <- df %>%
    filter(organism_category == organism) %>%
    mutate(
      successes = n_hospitalizations,
      failures = .data[[denominator_col]] - successes,
      exposure_value = .data[[exposure]]
    ) %>%
    filter(
      !is.na(exposure_value),
      !is.na(successes),
      !is.na(failures),
      .data[[denominator_col]] > 0,
      failures >= 0
    )

  if (nrow(dat) < 10 || sum(dat$successes, na.rm = TRUE) < 10 || sd(dat$exposure_value, na.rm = TRUE) == 0) {
    return(NULL)
  }

  exposure_iqr <- IQR(dat$exposure_value, na.rm = TRUE)
  if (!is.finite(exposure_iqr) || exposure_iqr <= 0) return(NULL)

  dat <- dat %>%
    mutate(exposure_iqr_scaled = exposure_value / exposure_iqr)

  fit <- tryCatch(
    glm(
      cbind(successes, failures) ~ exposure_iqr_scaled + factor(year),
      family = quasibinomial(),
      data = dat
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  coef_tab <- summary(fit)$coefficients
  if (!"exposure_iqr_scaled" %in% rownames(coef_tab)) return(NULL)

  beta <- coef_tab["exposure_iqr_scaled", "Estimate"]
  se <- coef_tab["exposure_iqr_scaled", "Std. Error"]
  p_value <- coef_tab["exposure_iqr_scaled", "Pr(>|t|)"]

  tibble(
    organism_category = organism,
    exposure = exposure,
    denominator = denominator_col,
    model_type = model_type,
    n_county_years = nrow(dat),
    total_detections = sum(dat$successes, na.rm = TRUE),
    total_denominator = sum(dat[[denominator_col]], na.rm = TRUE),
    exposure_iqr = exposure_iqr,
    odds_ratio_per_iqr = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    log_or = beta,
    se = se,
    p_value = p_value
  )
}

out_dir <- file.path("output", "final")

site_county_year_path <- latest_file(file.path(out_dir, "microbe_site_county_year_*.csv"))
organism_category_path <- latest_file(file.path(out_dir, "microbe_organism_category_*.csv"))

message("Using site-count county-year file: ", site_county_year_path)
message("Using organism-category file: ", organism_category_path)

site_county_year <- readr::read_csv(site_county_year_path, show_col_types = FALSE) %>%
  mutate(
    county_fips = str_pad(as.character(county_fips), 5, side = "left", pad = "0"),
    year = as.integer(year)
  )

organism_category <- readr::read_csv(organism_category_path, show_col_types = FALSE) %>%
  mutate(
    county_fips = str_pad(as.character(county_fips), 5, side = "left", pad = "0"),
    year = as.integer(year),
    organism_category = as.character(organism_category)
  )

analysis_dat <- tidyr::crossing(
  site_county_year %>%
    select(site_name, county_fips, year, n_icu, n_with_resp_culture,
           n_positive_resp_culture, pm25_mean, no2_mean),
  organism_category = sort(unique(organism_category$organism_category))
) %>%
  left_join(
    organism_category %>%
      select(site_name, county_fips, year, organism_category, n_hospitalizations, n_organism_isolates),
    by = c("site_name", "county_fips", "year", "organism_category")
  ) %>%
  mutate(
    n_hospitalizations = coalesce(n_hospitalizations, 0),
    n_organism_isolates = coalesce(n_organism_isolates, 0)
  )

eligible_organisms <- analysis_dat %>%
  group_by(organism_category) %>%
  summarise(total_detections = sum(n_hospitalizations, na.rm = TRUE), .groups = "drop") %>%
  filter(total_detections >= 25) %>%
  pull(organism_category)

model_grid <- tidyr::expand_grid(
  organism_category = eligible_organisms,
  exposure = c("pm25_mean", "no2_mean"),
  denominator_col = c("n_with_resp_culture", "n_positive_resp_culture")
) %>%
  mutate(
    model_type = case_when(
      denominator_col == "n_with_resp_culture" ~ "risk_among_respiratory_cultured",
      denominator_col == "n_positive_resp_culture" ~ "composition_among_positive_cultures",
      TRUE ~ denominator_col
    )
  )

model_results <- purrr::pmap_dfr(
  model_grid,
  ~ fit_one_model(
    df = analysis_dat,
    organism = ..1,
    exposure = ..2,
    denominator_col = ..3,
    model_type = ..4
  )
) %>%
  group_by(exposure, denominator, model_type) %>%
  mutate(fdr_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(exposure, denominator, fdr_p_value, desc(odds_ratio_per_iqr))

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- file.path(out_dir, glue("pollution_microbe_risk_models_{stamp}.csv"))
readr::write_csv(model_results, out_path)

message("Wrote model results: ", out_path)
message("")
message("Increased organism risk among respiratory-cultured hospitalizations, FDR < 0.10:")
model_results %>%
  filter(
    denominator == "n_with_resp_culture",
    odds_ratio_per_iqr > 1,
    fdr_p_value < 0.10
  ) %>%
  arrange(fdr_p_value, desc(odds_ratio_per_iqr)) %>%
  select(organism_category, exposure, total_detections, exposure_iqr,
         odds_ratio_per_iqr, ci_low, ci_high, p_value, fdr_p_value) %>%
  print(n = 50)

message("")
message("Top positive associations among respiratory-cultured hospitalizations by nominal p-value:")
model_results %>%
  filter(
    denominator == "n_with_resp_culture",
    odds_ratio_per_iqr > 1
  ) %>%
  arrange(p_value) %>%
  select(organism_category, exposure, total_detections, exposure_iqr,
         odds_ratio_per_iqr, ci_low, ci_high, p_value, fdr_p_value) %>%
  head(20) %>%
  print(n = 20)
