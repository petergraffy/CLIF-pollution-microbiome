#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(arrow);library(dplyr);library(jsonlite)})
source('utils/mwas.R');source('utils/mwas_power.R');source('utils/clif_io.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
b <- readRDS(file.path(run_dir,'private','cohort.rds'))
m <- readRDS(file.path(run_dir,'private','matched_exposures.rds'))
det <- b$detections[hours_from_icu<=48]
orgs <- sort(unique(det$organism_category)); family_size <- length(orgs)*2L
rows <- list();i <- 0L
for(org in orgs) for(ex in c('pm25','o3')) {
  ids <- unique(det[organism_category==org,hospitalization_id])
  info <- mwas_design_information(m[stratum %in% ids],paste0(ex,'_lag1_7'),if(ex=='pm25')5 else 10)
  for(a in c(.05,.05/family_size)) for(target in c(1.05,1.1,1.2,1.5,2)) {
    need <- mwas_required_information(log(target),a)
    i <- i+1L
    rows[[i]] <- as.data.table(c(list(organism=org,pollutant=ex,
      exposure_unit=if(ex=='pm25')'5 ug/m3' else '10 ppb',n_eligible_events=length(ids)),info,
      list(alpha=a,multiplicity=if(a==.05)'nominal' else 'Bonferroni_all_214_or_actual_family',
      family_size=family_size,target_or=target,target_power=.8,
      approximate_power=mwas_normal_power(log(target),info$information,a),
      approximate_events_required=if(is.finite(info$information)&&info$information>1e-10)
        ceiling(need/(info$information/info$n_events)) else NA_real_,
      approximate_mde_or_80=if(is.finite(info$information)&&info$information>1e-10)
        exp(sqrt(mwas_required_information(1,a)/info$information)) else NA_real_)))
  }
}
fwrite(rbindlist(rows),file.path(run_dir,'power_design_scenarios_48h.csv'))
# Original source labels: aggregate only, no identifiers written to the audit.
ids <- b$cohort$hospitalization_id
raw <- as.data.table(open_dataset(find_table_path('microbiology_culture')) %>%
  select(hospitalization_id,order_dttm,collect_dttm,fluid_category,fluid_name,
         method_category,method_name,lab_loinc_code) %>%
  filter(hospitalization_id %in% ids) %>% collect())
raw[,hospitalization_id:=as.character(hospitalization_id)]
for(v in c('order_dttm','collect_dttm'))set(raw,j=v,value=mwas_ts(raw[[v]]))
for(v in c('fluid_category','method_category','method_name'))set(raw,j=v,value=mwas_clean(raw[[v]]))
raw[,episode_id:=mwas_episode_key(.SD)]
raw <- merge(raw,b$episodes[hours_from_icu<=48,.(episode_id)],by='episode_id')
audit <- unique(raw[,.(episode_id,hospitalization_id,fluid_category,fluid_name)])[
  ,.(n_episodes=uniqueN(episode_id),n_hospitalizations=uniqueN(hospitalization_id)),
  by=.(fluid_category,fluid_name)]
setorder(audit,-n_episodes)
fwrite(audit,file.path(run_dir,'specimen_source_audit_48h.csv'))
write_json(list(method='Null conditional likelihood information, weather/holiday residualized; local normal power approximation',
  culture_window_hours=48,family_size=family_size,target_power=.8,target_or=c(1.05,1.1,1.2,1.5,2),
  limitations=c('Not observed-effect post hoc power','Independent matched sets assumed; patient and shared calendar exposure dependence not simulated',
  'Uniform case probabilities assume zero nuisance coefficients; non-null weather effects not modeled',
  'Required event count extrapolates the existing organism-specific exposure distribution',
  'Bonferroni is conservative family-wise control, not an estimate of BH discovery power',
  'Sparse designs and large effects require simulation; no eligibility changes made'),
  specimens='Existing broad four-category definition audited, not changed'),
  file.path(run_dir,'power_audit_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
message('Wrote power scenarios for ',length(orgs),' organisms and original specimen source audit')
