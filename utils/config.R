if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required to load config/config.json.")
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

load_config <- function(config_path = Sys.getenv("CLIF_CONFIG_PATH", unset = "config/config.json"),
                        required = TRUE) {
  if (file.exists(config_path)) {
    message("Loaded configuration from ", config_path)
    site_config <- jsonlite::fromJSON(config_path)
    # Shared protocol defaults are supplied centrally, not entered by each site.
    defaults <- list(site_timezone = "America/Chicago",
                     mwas_exposure_cache = "data/public/exposures",
                     mwas_acs_dir = "data/public/acs/2017", derive_sofa = FALSE)
    return(utils::modifyList(defaults, site_config))
  }

  if (required) {
    stop(
      "Configuration file not found at ", config_path,
      ". Copy config/config_template.json to config/config.json and update local paths, ",
      "or set CLIF_CONFIG_PATH."
    )
  }

  list()
}

config_value <- function(config, fields, env = NULL, default = NULL, required = FALSE) {
  if (!is.null(env)) {
    env_value <- Sys.getenv(env, unset = NA_character_)
    if (!is.na(env_value) && nzchar(env_value)) return(env_value)
  }

  for (field in fields) {
    value <- config[[field]]
    if (!is.null(value) && length(value) > 0 && !all(is.na(value)) && !identical(value, "")) {
      return(value)
    }
  }

  if (required && is.null(default)) {
    stop("Missing required configuration field: ", paste(fields, collapse = " or "))
  }

  default
}

config <- load_config(required = TRUE)
