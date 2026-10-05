# Local site configuration

Copy `config_template.json` to `config.json`, which is ignored by Git. Set a unique site name, the CLIF 2.1 release folder, `file_type="parquet"`, and the hospital's IANA timezone. Actual calendar dates and correctly interpreted timestamp offsets are required.

`mwas_exposure_cache` defaults to `data/public/exposures`; `mwas_acs_dir` defaults to `data/public/acs/2017`. These are bundled national public inputs. Keep those defaults for the shared buddy protocol. Paths are relative to the repository root. `MWAS_EXPOSURE_CACHE` and `MWAS_ACS_DIR` override these fields for deliberate alternate inputs. `CLIF_CONFIG_PATH` selects an alternate local configuration file.

The six core tables are patient, hospitalization, adt, microbiology_culture, medication_admin_intermittent and hospital_diagnosis. Run `Rscript code/33_site_preflight.R` to check required fields. Ethnicity and admission type are optional descriptive fields. Modified SOFA additionally requires labs, vitals, respiratory support, assessments and continuous medication administrations; enable with `MWAS_SKIP_SEVERITY=0`.

Ordinary site runs make no environmental network requests and require no Census key. Optional public ACS rebuilds can use the summary-file route without a key; the API alternative reads `CENSUS_API_KEY` only from the local environment. Never commit credentials, completed configurations, source tables or clinical linkage files.
