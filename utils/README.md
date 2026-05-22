## Utils directory

This directory contains shared helper code used by the project scripts.

Current utilities:

1. `config.R`: loads `config/config.json` for R scripts.
2. `config.py`: parallel Python config loader retained from the CLIF project template.

The active workflow is currently R-based. Scripts source `utils/config.R` to access `site_name`, `tables_path`, `file_type`, `repo`, and `exposome_path`.
