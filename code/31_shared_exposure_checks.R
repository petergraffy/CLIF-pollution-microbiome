#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
source('utils/mwas_federated.R');source('utils/mwas_poisson.R');source('utils/mwas_targets.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
root <- file.path(run_dir,'federated');input <- fread(file.path(root,'site_estimates.csv'))
b <- readRDS(file.path(run_dir,'private','cohort.rds'));m <- readRDS(file.path(run_dir,'private','matched_exposures.rds'))
groups <- readRDS(file.path(run_dir,'private','admission_diagnosis_groups.rds'))
m <- merge(m,groups,by.x='stratum',by.y='hospitalization_id',all.x=TRUE,sort=FALSE)
common <- setNames(lapply(mwas_adjustment_labels,function(adj)
 setNames(lapply(c('pm25','o3'),function(ex)mwas_common_windows(m,ex,weather_adjusted=adj=='weather_adjusted')),c('pm25','o3'))),mwas_adjustment_labels)
targets <- input[term=='pollution' & !startsWith(analysis,'ses:') & !startsWith(analysis,'modifier:')]
rows <- vector('list',nrow(targets));checks <- vector('list',nrow(targets))
for(i in seq_len(nrow(targets))) {
 if(i%%500L==0L)message('Shared-exposure checks: ',i,'/',nrow(targets))
 r <- targets[i];d <- common[[r$model_adjustment]][[r$pollutant]][stratum %in% mwas_target_ids(b,r$organism,r$analysis)]
 if(r$analysis %in% c('pneumonia_aspiration','obstructive_airway','other_respiratory','nonrespiratory'))d <- d[diagnosis_group==r$analysis]
 f <- list(status=paste0('not_attempted_',r$status))
 if(r$status=='ok') {
  f <- tryCatch(mwas_count_fit(mwas_count_design(d,paste0(r$pollutant,'_',r$exposure_window),
     if(r$pollutant=='pm25')5 else 10,weather_adjusted=r$model_adjustment=='weather_adjusted')),error=function(e)list(status='fit_error'))
 }
 copies <- rbind(copy(r),copy(r));copies[,inference_method:=c('conditional_poisson_quasi','conditional_poisson_hac28')]
 copies[,`:=`(status=f$status,log_or=NA_real_,se=NA_real_,null_information=NA_real_)]
 if(f$status=='ok')copies[,`:=`(log_or=f$log_or,se=unname(f$se[c('quasi','hac')]))]
 rows[[i]] <- copies
 checks[[i]] <- data.table(organism=r$organism,pollutant=r$pollutant,analysis=r$analysis,exposure_window=r$exposure_window,model_adjustment=r$model_adjustment,
  status=f$status,clogit_status=r$status,n_events=r$n_events,
  point_difference=if(f$status=='ok')f$log_or-r$log_or else NA_real_,
  dispersion=if(f$status=='ok')f$dispersion else NA_real_,
  n_calendar_days=if(f$status=='ok')f$n_calendar_days else NA_integer_,
  quasi_se=if(f$status=='ok')f$se['quasi'] else NA_real_,hac_se=if(f$status=='ok')f$se['hac'] else NA_real_,
  patient_se=r$se)
}
fwrite(rbindlist(rows),file.path(root,'shared_exposure_estimates.csv'))
fwrite(rbindlist(checks),file.path(root,'shared_exposure_diagnostics.csv'))
message('Conditional Poisson variance checks complete: ',nrow(targets),' attempted models')
