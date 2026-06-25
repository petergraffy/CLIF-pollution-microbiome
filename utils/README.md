## Utils directory

This directory contains shared helper code used by the project scripts.

Current utilities:

1. `config.R`: loads `config/config.json` for R scripts.
2. `clif_io.R`: shared CLIF table discovery and ZCTA exposure file helpers.
3. `config.py`: parallel Python config loader retained from the CLIF project template.

The active workflow is currently R-based. SHRF/ZCTA scripts source `utils/clif_io.R` to access `site_name`, `tables_path`, `file_type`, `repo`, and `zcta_exposure_dir`.
