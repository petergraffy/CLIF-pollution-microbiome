#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(arrow);library(dplyr);library(jsonlite)})
source('utils/clif_io.R');source('utils/mwas.R');source('utils/mwas_federated.R');source('utils/mwas_targets.R');source('utils/mwas_characteristics.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
acs_dir <- config_value(config,'mwas_acs_dir',env='MWAS_ACS_DIR',default='data/public/acs/2017')
exposure_dir <- config_value(config,'mwas_exposure_cache',env='MWAS_EXPOSURE_CACHE',default='data/public/exposures')
exposure_manifest_path <- file.path(exposure_dir,'manifest.json')
out_dir <- file.path(run_dir,'federated');dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
b <- readRDS(file.path(run_dir,'private','cohort.rds'))
m <- readRDS(file.path(run_dir,'private','matched_exposures.rds'))
characteristics_path <- file.path(run_dir,'private','site_characteristics.rds')
if(!file.exists(characteristics_path))stop('Run code/37_site_characteristics.R before the federated models')
characteristics <- readRDS(characteristics_path)
modifiers <- mwas_modifier_registry();fwrite(modifiers,file.path(out_dir,'modifier_registry.csv'))
modifier_columns <- unique(modifiers$column)
for(v in setdiff(modifier_columns,names(characteristics)))characteristics[,(v):=NA_real_]
linked <- characteristics[,c('hospitalization_id',modifier_columns),with=FALSE]
m <- merge(m,linked,by.x='stratum',by.y='hospitalization_id',all.x=TRUE,sort=FALSE)
windows <- c(3L,7L,14L,28L)
# Primary sets do not require weather; sensitivity sets require complete weather.
common <- setNames(lapply(mwas_adjustment_labels,function(adj)
 setNames(lapply(c('pm25','o3'),function(ex)mwas_common_windows(m,ex,windows,adj=='weather_adjusted')),c('pm25','o3'))),mwas_adjustment_labels)
manifest <- read_json(file.path(run_dir,'manifest.json'),simplifyVector=TRUE)
acs_path <- file.path(acs_dir,'zcta_ses.csv')
have_acs <- file.exists(acs_path)
registry <- data.table(indicator=c('poverty_pct','no_high_school_pct','unemployment_pct','median_household_income','crowding_pct'),
 center=c(20,20,5,50000,2),scale=c(10,10,5,NA,5),transform=c('linear','linear','linear','log2','linear'),
 family=c('ses_primary','ses_secondary','ses_secondary','ses_secondary','ses_secondary'))
fwrite(registry,file.path(out_dir,'ses_registry.csv'))
acs_manifest <- if(have_acs)read_json(file.path(acs_dir,'manifest.json'),simplifyVector=TRUE) else NULL
if(have_acs) {
 acs <- fread(acs_path,colClasses=list(character='zcta'));stopifnot(!anyDuplicated(acs$zcta),uniqueN(acs$acs_year)==1)
 m <- merge(m,acs,by.x='zip',by.y='zcta',all.x=TRUE,sort=FALSE)
} else message('ACS cache missing: overall and diagnosis models will run; SES models explicitly unavailable')
ids <- unique(b$episodes[hours_from_icu<=48,hospitalization_id])
capture <- rbindlist(lapply(seq_len(nrow(modifiers)),function(j) {
 reg <- modifiers[j];values <- mwas_modifier_values(characteristics[hospitalization_id %in% ids],reg)
 data.table(indicator=reg$indicator,n_admissions=length(values),n_observed=sum(is.finite(values)),n_excluded=sum(!is.finite(values)),
 n_reference=if(reg$type=='contrast')sum(values==0,na.rm=TRUE) else NA_integer_,
 n_comparison=if(reg$type=='contrast')sum(values==1,na.rm=TRUE) else NA_integer_)
}))
fwrite(capture,file.path(out_dir,'modifier_capture.csv'))
window_qc <- rbindlist(lapply(mwas_adjustment_labels,function(adj)
 rbindlist(lapply(c('pm25','o3'),function(ex)rbindlist(lapply(windows,function(w) {
 single <- mwas_common_windows(m[stratum %in% ids],ex,w,adj=='weather_adjusted')
 shared <- common[[adj]][[ex]][stratum %in% ids]
 data.table(model_adjustment=adj,pollutant=ex,window_days=w,n_window_available_events=uniqueN(single$stratum),
   n_common_events=uniqueN(shared$stratum),n_common_referents=sum(shared$case==0L))
}))))))
fwrite(window_qc,file.path(out_dir,'exposure_window_qc.csv'))
dx <- as.data.table(open_dataset(find_table_path('hospital_diagnosis')) %>%
 select(hospitalization_id,diagnosis_code,diagnosis_code_format,diagnosis_primary,poa_present) %>%
 filter(hospitalization_id %in% ids) %>% collect())
groups <- mwas_admission_diagnoses(dx,ids)
saveRDS(groups,file.path(run_dir,'private','admission_diagnosis_groups.rds'))
fwrite(groups[,.(n_admissions=.N),by=diagnosis_group],file.path(out_dir,'diagnosis_capture.csv'))
m <- merge(m,groups,by.x='stratum',by.y='hospitalization_id',all.x=TRUE,sort=FALSE)
# This QC table contains aggregate linkage counts only, not ZIP codes.
if(have_acs)fwrite(rbindlist(lapply(registry$indicator,function(v)
 m[case==1 & stratum %in% ids,.(indicator=v,n_admissions=.N,n_observed=sum(is.finite(get(v))),n_missing=sum(!is.finite(get(v))))])),
 file.path(out_dir,'ses_linkage_qc.csv'))
det <- b$detections[hours_from_icu<=72]
orgs <- c(sort(unique(det$organism_category)),'__early_icu_admission','__any_respiratory_culture','__any_named_organism')
protocol <- list(version='acute_mwas_federated_v6_weather_sensitivity',culture_window=48,icu_entry_hours=24,
 source_categories=c('respiratory_tract','respiratory_tract_lower','nasopharynx_upperairway','oropharynx_tongue_oralcavity'),
 method='culture',event='local hospital admission date',lag='mean complete days 1-N',
 exposure_windows=windows,primary_exposure_window=7,
 comparison_population='identical exposure-complete/referent rows across all four windows within pollutant and adjustment specification; sensitivity additionally requires weather; nonzero exposure variation in every window',
 matching='same weekday, month, year, residential ZIP proxy',primary_adjustment='holiday only; no weather terms',
 weather_sensitivity='ns seven-day temperature df3 + ns seven-day RH df3 + holiday; all analysis families and windows',
 model_adjustments=mwas_adjustment_labels,
 inference='Primary Efron conditional logistic patient covariance; conditional Poisson quasi and calendar HAC28 checks',
 clinical_sensitivities=c('culture24','culture72','no_documented_antibiotics','pulmonary_sources','upper_airway_sources'),
 selection_companions=c('early_icu_admission','any_respiratory_culture','any_named_organism'),
 age_restriction='none',
 effect_modifiers=as.data.frame(modifiers),modifier_window=7,
 modifier_method='separate pollution interactions; no modifier main effects; fixed race category contrasts versus white; unknown/missing excluded; complete SOFA only; 24h sensitivity',
 model_package_versions=as.list(setNames(vapply(c('survival','data.table','lubridate'),function(p)as.character(packageVersion(p)),character(1)),c('survival','data.table','lubridate'))),
 start_date=manifest$start_date,end_date=manifest$end_date,
 public_exposure_manifest_md5=if(file.exists(exposure_manifest_path))unname(tools::md5sum(exposure_manifest_path)) else 'alternate_unmanifested_fixture',
 acs_year=if(have_acs)unique(acs$acs_year) else NA_integer_,
 acs_indicator_md5=if(have_acs)unname(tools::md5sum(acs_path)) else 'unavailable',
 antibacterial_mapping_md5=manifest$antibiotic_mapping_md5,
 code_md5=as.list(tools::md5sum(c('utils/mwas.R','utils/mwas_federated.R','utils/mwas_targets.R',
  'utils/mwas_poisson.R','utils/mwas_characteristics.R','code/37_site_characteristics.R','code/27_run_federated_mwas.R','code/31_shared_exposure_checks.R',
  'code/21_prepare_acute_mwas.R','code/22_cache_mwas_exposures.py','code/23_run_acute_mwas.R'))),
 diagnosis='primary and POA==1; grouped ICD9/10; missing/ambiguous not nonrespiratory',ses=as.data.frame(registry))
write_json(protocol,file.path(out_dir,'protocol.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
protocol_id <- unname(tools::md5sum(file.path(out_dir,'protocol.json')))
rows <- list();i <- 0L
common <- setNames(lapply(mwas_adjustment_labels,function(adj)
 setNames(lapply(c('pm25','o3'),function(ex)m[stratum %in% common[[adj]][[ex]]$stratum & row_id %in% common[[adj]][[ex]]$row_id]),c('pm25','o3'))),mwas_adjustment_labels)
analyses <- c('overall','pneumonia_aspiration','obstructive_airway','other_respiratory','nonrespiratory',paste0('ses:',registry$indicator),
 paste0('clinical:',c('culture24','culture72','no_documented_antibiotics','pulmonary_sources','upper_airway_sources')),modifiers$analysis)
for(org in orgs) {
 selected_analyses <- if(startsWith(org,'__'))'selection_companion' else analyses
 for(adj in mwas_adjustment_labels)for(analysis in selected_analyses)for(ex in c('pm25','o3'))for(w in if(startsWith(analysis,'clinical:') || startsWith(analysis,'modifier:'))7L else windows) {
  org_ids <- mwas_target_ids(b,org,analysis)
  d <- copy(common[[adj]][[ex]][stratum %in% org_ids]); interaction <- startsWith(analysis,'ses:')
  modifier <- startsWith(analysis,'modifier:')
  modifier_index <- if(modifier)match(analysis,modifiers$analysis) else NA_integer_
  modifier_reg <- if(modifier)modifiers[modifier_index] else NULL
  indicator <- if(modifier)modifier_reg$indicator else if(interaction)sub('^ses:','',analysis) else 'none'
  family <- if(analysis=='overall')'primary' else if(!interaction)'diagnosis_secondary' else NA_character_
  if(startsWith(analysis,'clinical:'))family <- 'clinical_sensitivity'
  if(analysis=='selection_companion')family <- 'selection_companion'
  if(modifier)family <- modifier_reg$family
  # Registry lookup without ambiguous data.table column/variable scoping.
  reg <- if(interaction)registry[match(sub('^ses:','',analysis),registry$indicator)] else NULL
  if(interaction)family <- reg$family
  if(w!=7L)family <- if(analysis=='overall')'exposure_duration_sensitivity' else paste0(family,'_duration_sensitivity')
  if(adj=='weather_adjusted')family <- paste0(family,'_weather_sensitivity')
  if(analysis %in% c('pneumonia_aspiration','obstructive_airway','other_respiratory','nonrespiratory'))d <- d[diagnosis_group==analysis]
  if(interaction && !have_acs) {
   f <- data.table(term=c('pollution','pollution_ses'),n_events=0L,n_patients=0L,n_referents=0L,
     design_rank=NA_integer_,status='acs_unavailable',log_or=NA_real_,se=NA_real_,cov_main_interaction=NA_real_,null_information=NA_real_)
  } else {
   if(interaction) {
    v <- d[[indicator]]
    ses_value <- if(reg$transform=='log2')ifelse(is.finite(v)&v>0,log2(v/reg$center),NA_real_) else (v-reg$center)/reg$scale
    d[,ses:=ses_value]
   }
   if(modifier)d[,ses:=mwas_modifier_values(d,modifier_reg)]
   f <- mwas_federated_fit(d,paste0(ex,'_lag1_',w),if(ex=='pm25')5 else 10,interaction || modifier,weather_adjusted=adj=='weather_adjusted')
  }
  if(modifier) {
   f[term=='pollution_ses',term:='pollution_modifier']
   f[status=='no_ses_variation',status:='no_modifier_variation']
  }
  f[,`:=`(modifier_reference=if(modifier)modifier_reg$reference else NA_character_,
    modifier_comparison=if(modifier)modifier_reg$comparison else NA_character_)]
  f[,`:=`(site=clif_site_name,protocol_id=protocol_id,organism=org,pollutant=ex,
    culture_window_hours=if(analysis=='clinical:culture24')24L else if(analysis=='clinical:culture72')72L else 48L,
    inference_method='patient_cluster',model_adjustment=adj,exposure_window=paste0('lag1_',w),window_days=w,
    exposure_unit=if(ex=='pm25')'5 ug/m3' else '10 ppb',
    analysis=analysis,family=family,indicator=indicator)]
  i <- i+1L;rows[[i]] <- f
 }
 message('Completed federated models: ',org)
}
res <- rbindlist(rows)
# Explicit export schema: no identifiers, ZIPs, dates, patient rows, or raw warning strings.
fwrite(res,file.path(out_dir,'site_estimates.csv'))
fwrite(res[,.(n_model_terms=.N,n_ok=sum(status=='ok')),by=.(analysis,exposure_window,model_adjustment,status)],file.path(out_dir,'model_diagnostics.csv'))
write_json(list(site=clif_site_name,protocol_id=protocol_id,acs_available=have_acs,
 acs_manifest=acs_manifest,named_organisms=sum(!startsWith(orgs,'__')),
 inference_limitations=c('Sparse estimates are exploratory even when they pass rank and covariance checks',
 'Patient-cluster covariance does not address all shared date/ZCTA exposure dependence',
 'Diagnosis strata define selected outcomes; not baseline confounder adjustment',
 'ACS sampling uncertainty not propagated; residential ZIP is a ZCTA proxy'),
 effect_modifier_limitations=c('Each modifier is fit separately; race is recorded social classification, not a biological mechanism',
 'Six-hour SOFA is early post-ICU organ dysfunction; 24h score is later and may be affected by exposure/treatment',
 'Complete-score and demographic missingness may select patients; modifier models currently use patient covariance only'),
 session_info=capture.output(sessionInfo())),file.path(out_dir,'site_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
message('Federated aggregate outputs: ',out_dir)
