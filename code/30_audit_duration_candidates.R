#!/usr/bin/env Rscript
# Local fit diagnostics for FDR-passing overall candidates; not a new screen.
suppressPackageStartupMessages({library(data.table);library(survival)})
source('utils/mwas_federated.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
root <- file.path(run_dir,'federated')
pool <- Sys.getenv('MWAS_POOL_DIR',file.path(root,'pool_ucmc'))
p <- fread(file.path(pool,'pooled_mwas.csv'))
candidates <- p[analysis=='overall' & inference_method=='patient_cluster' & q_value<.05 & k_sites==1]
if(!nrow(candidates)) {
 stale <- file.path(root,'duration_candidate_diagnostics.csv')
 if(file.exists(stale))unlink(stale)
 message('No single-site overall candidates to audit');quit(save='no')
}
b <- readRDS(file.path(run_dir,'private','cohort.rds'))
m <- readRDS(file.path(run_dir,'private','matched_exposures.rds'))
rows <- lapply(seq_len(nrow(candidates)),function(i) {
 r <- candidates[i]
 ids <- unique(b$detections[hours_from_icu<=48 & organism_category==r$organism,hospitalization_id])
 weather_adjusted <- r$model_adjustment=='weather_adjusted'
 d <- mwas_common_windows(m[stratum %in% ids],r$pollutant,weather_adjusted=weather_adjusted)
 exposure <- paste0(r$pollutant,'_',r$exposure_window)
 unit <- if(r$pollutant=='pm25')5 else 10
 d[,x:=get(exposure)/unit]
 f <- survival::clogit(mwas_model_formula(weather_adjusted=weather_adjusted),data=d,method='efron',control=coxph.control(iter.max=50))
 idx <- match('x',names(coef(f)));model_se <- sqrt(f$naive.var[idx,idx])
 loo <- rbindlist(lapply(unique(d$patient_id),function(h)
   mwas_federated_fit(d[patient_id!=h],exposure,unit,weather_adjusted=weather_adjusted)))
 data.table(organism=r$organism,pollutant=r$pollutant,exposure_window=r$exposure_window,model_adjustment=r$model_adjustment,
  n_events=r$n_events,log_or=r$log_or,robust_se=r$se,model_based_se=model_se,
  model_based_p=2*pnorm(-abs(r$log_or/model_se)),q_value=r$q_value,
  leave_one_patient_out_attempted=nrow(loo),leave_one_patient_out_ok=sum(loo$status=='ok'),
  leave_one_patient_out_failed=sum(loo$status!='ok'),
  review='Exploratory candidate; inspect effect size, uncertainty and refit failures before interpretation')
})
fwrite(rbindlist(rows),file.path(root,'duration_candidate_diagnostics.csv'))
message('Audited ',length(rows),' local overall candidates')
