#!/usr/bin/env Rscript
# Aggregate support for every interaction plus sparse-candidate stability checks.
suppressPackageStartupMessages({library(data.table);library(survival)})
source('utils/mwas.R');source('utils/mwas_federated.R');source('utils/mwas_targets.R');source('utils/mwas_characteristics.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
root <- file.path(run_dir,'federated');pool <- Sys.getenv('MWAS_POOL_DIR',file.path(root,'pool_ucmc'))
input <- fread(file.path(root,'site_estimates.csv'))
b <- readRDS(file.path(run_dir,'private','cohort.rds'));m <- readRDS(file.path(run_dir,'private','matched_exposures.rds'))
h <- readRDS(file.path(run_dir,'private','site_characteristics.rds'));registry <- mwas_modifier_registry()
columns <- unique(registry$column);m <- merge(m,h[,c('hospitalization_id',columns),with=FALSE],by.x='stratum',by.y='hospitalization_id',all.x=TRUE,sort=FALSE)
common <- setNames(lapply(mwas_adjustment_labels,function(adj)
 setNames(lapply(c('pm25','o3'),function(ex)mwas_common_windows(m,ex,weather_adjusted=adj=='weather_adjusted')),c('pm25','o3'))),mwas_adjustment_labels)
cases <- lapply(common,function(spec)lapply(spec,function(d)d[case==1L]));targets <- input[startsWith(analysis,'modifier:') & term=='pollution_modifier']
support <- rbindlist(lapply(seq_len(nrow(targets)),function(i) {
 r <- targets[i];j <- match(r$analysis,registry$analysis);reg <- registry[j]
 d <- cases[[r$model_adjustment]][[r$pollutant]][stratum %in% mwas_target_ids(b,r$organism)]
 v <- mwas_modifier_values(d,reg);ok <- is.finite(v)
 data.table(organism=r$organism,pollutant=r$pollutant,analysis=r$analysis,exposure_window=r$exposure_window,model_adjustment=r$model_adjustment,
  n_before_modifier_exclusions=nrow(d),n_events_with_modifier=sum(ok),n_excluded=sum(!ok),
  n_reference=if(reg$type=='contrast')sum(v==0,na.rm=TRUE) else NA_integer_,
  n_comparison=if(reg$type=='contrast')sum(v==1,na.rm=TRUE) else NA_integer_,
  n_patients_with_modifier=uniqueN(d$patient_id[ok]),
  n_reference_patients=if(reg$type=='contrast')uniqueN(d$patient_id[!is.na(v)&v==0]) else NA_integer_,
  n_comparison_patients=if(reg$type=='contrast')uniqueN(d$patient_id[!is.na(v)&v==1]) else NA_integer_,
  n_distinct_modifier_values=uniqueN(v[ok]),
  modifier_min=if(any(ok))min(v[ok]) else NA_real_,modifier_max=if(any(ok))max(v[ok]) else NA_real_,
  site=r$site,protocol_id=r$protocol_id)
}))
fwrite(support,file.path(root,'modifier_case_support.csv'))
p <- fread(file.path(pool,'pooled_mwas.csv'))
candidates <- p[startsWith(analysis,'modifier:') & term=='pollution_modifier' & inference_method=='patient_cluster' & q_value<.05 & k_sites==1]
path <- file.path(root,'modifier_candidate_diagnostics.csv')
if(!nrow(candidates)) {if(file.exists(path))unlink(path);message('No local modifier candidates; all model support exported');quit(save='no')}
rows <- lapply(seq_len(nrow(candidates)),function(i) {
 r <- candidates[i];j <- match(r$analysis,registry$analysis);reg <- registry[j]
 d <- copy(common[[r$model_adjustment]][[r$pollutant]][stratum %in% mwas_target_ids(b,r$organism)])
 d[,ses:=mwas_modifier_values(d,reg)];d <- d[is.finite(ses)]
 exposure <- paste0(r$pollutant,'_',r$exposure_window);unit <- if(r$pollutant=='pm25')5 else 10
 d[,`:=`(x=get(exposure)/unit,xs=get(exposure)/unit*ses)]
 warnings <- character()
 fit <- tryCatch(withCallingHandlers(survival::clogit(mwas_model_formula(TRUE,r$model_adjustment=='weather_adjusted'),data=d,method='efron',control=coxph.control(iter.max=50)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)NULL)
 model_se <- if(!is.null(fit))sqrt(fit$naive.var[match('xs',names(coef(fit))),match('xs',names(coef(fit)))]) else NA_real_
 loo <- rbindlist(lapply(unique(d$patient_id),function(id) {
  f <- mwas_federated_fit(d[patient_id!=id],exposure,unit,TRUE,weather_adjusted=r$model_adjustment=='weather_adjusted');f[term=='pollution_ses']
 }))
 row <- data.table(organism=r$organism,pollutant=r$pollutant,analysis=r$analysis,exposure_window=r$exposure_window,model_adjustment=r$model_adjustment,
  n_events=r$n_events,log_ratio_of_ORs=r$log_or,patient_cluster_se=r$se,model_based_se=model_se,
  model_based_p=if(is.finite(model_se))2*pnorm(-abs(r$log_or/model_se)) else NA_real_,q_value=r$q_value,
  leave_one_patient_out_attempted=nrow(loo),leave_one_patient_out_ok=sum(loo$status=='ok'),leave_one_patient_out_failed=sum(loo$status!='ok'),
  loo_min_log_ratio=if(any(loo$status=='ok'))min(loo[status=='ok',log_or]) else NA_real_,
  loo_max_log_ratio=if(any(loo$status=='ok'))max(loo[status=='ok',log_or]) else NA_real_,
  naive_fit_warned=length(warnings)>0,
  review='Exploratory only: inspect group support, model-based uncertainty and refit failures; baseline simulations do not validate modifier interactions')
 merge(row,support[,!c('site','protocol_id'),with=FALSE],by=c('organism','pollutant','analysis','exposure_window','model_adjustment'),all.x=TRUE)
})
fwrite(rbindlist(rows),path);message('Audited ',nrow(candidates),' local modifier candidates')
