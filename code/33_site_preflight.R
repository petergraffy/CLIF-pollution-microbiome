#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(jsonlite);library(arrow)})
source('utils/clif_io.R')
required <- list(
 patient=c('patient_id','sex_category','race_category'),
 hospitalization=c('patient_id','hospitalization_id','admission_dttm','discharge_dttm','age_at_admission','zipcode_five_digit','discharge_category'),
 adt=c('hospitalization_id','in_dttm','out_dttm','location_category'),
 microbiology_culture=c('hospitalization_id','organism_id','order_dttm','collect_dttm','result_dttm','fluid_category','fluid_name',
  'organism_name','method_category','method_name','organism_category','organism_group','lab_loinc_code'),
 medication_admin_intermittent=c('hospitalization_id','admin_dttm','med_category','med_route_category','mar_action_group'),
 hospital_diagnosis=c('hospitalization_id','diagnosis_code','diagnosis_code_format','diagnosis_primary','poa_present'))
tz <- config_value(config,'site_timezone',env='MWAS_TIMEZONE',required=TRUE)
if(!tz %in% OlsonNames())stop('Invalid site_timezone; supply the hospital IANA timezone')
if(!nzchar(clif_site_name)||clif_site_name=='Your_Site_Name')stop('Set a unique real site_name')
for(table in names(required)) {
 path <- find_table_path(table);columns <- names(arrow::open_dataset(path))
 missing <- setdiff(required[[table]],columns)
 if(length(missing))stop(table,' missing columns: ',paste(missing,collapse=', '))
 message(table,': schema OK')
}
if(!file.exists('resources/mwas/clif_intermittent_med_categories.csv'))stop('Bundled public medication vocabulary is missing')
for(pkg in c('data.table','survival','jsonlite','dplyr','lubridate','comorbidity','ggplot2'))
 if(!requireNamespace(pkg,quietly=TRUE))stop('Missing R dependency: ',pkg)
message('Core site preflight passed. Verify unshifted calendar dates and correct timezone with your data steward before analysis.')
