# ================================================================================================
# CLIF Pollution-Microbiome | Hierarchical Sensitivity Summary
# Purpose:
#   Summarize main and sensitivity hierarchical model outputs.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(glue)
})

parse_label <- function(path) {
  base <- basename(path)
  label <- str_match(base, "hierarchical_pollution_microbe_phenotype_models_UCMC_(.*)_\\d{8}_\\d{6}\\.csv")[, 2]
  ifelse(is.na(label), "legacy_all_counties", label)
}

files <- Sys.glob("output/final/hierarchical_pollution_microbe_phenotype_models_UCMC*.csv")
if (length(files) == 0) stop("No hierarchical model files found.")

all_results <- purrr::map_dfr(files, function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      analysis_label = parse_label(path),
      source_file = basename(path)
    )
})

summary <- all_results %>%
  group_by(analysis_label, source_file) %>%
  summarise(
    n_rows = n(),
    n_organisms = n_distinct(organism_category),
    n_phenotype_present_fdr_0_10_positive = sum(
      phenotype_level == "present" &
        odds_ratio_per_iqr > 1 &
        fdr_p_value < 0.10,
      na.rm = TRUE
    ),
    n_nominal_positive = sum(
      phenotype_level == "present" &
        odds_ratio_per_iqr > 1 &
        p_value < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(analysis_label, source_file)

signals <- all_results %>%
  filter(
    phenotype_level == "present",
    odds_ratio_per_iqr > 1,
    fdr_p_value < 0.10
  ) %>%
  arrange(analysis_label, fdr_p_value) %>%
  select(
    analysis_label, organism_category, exposure, phenotype,
    n_detected, n_phenotype_present, odds_ratio_per_iqr,
    ci_low, ci_high, p_value, fdr_p_value,
    interaction_or, interaction_p_value, source_file
  )

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
summary_path <- file.path("output", "final", glue("hierarchical_sensitivity_summary_{stamp}.csv"))
signals_path <- file.path("output", "final", glue("hierarchical_sensitivity_fdr_signals_{stamp}.csv"))

readr::write_csv(summary, summary_path)
readr::write_csv(signals, signals_path)

message("Wrote sensitivity summary: ", summary_path)
message("Wrote FDR signal table: ", signals_path)
message("")
print(summary, n = nrow(summary))
message("")
print(signals, n = min(nrow(signals), 100))
