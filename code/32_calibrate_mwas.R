#!/usr/bin/env Rscript
# Conditional simulation using locally observed exposure series and matched sets.
# Expanded patient scores and count/HAC scores share the SAME likelihood fit.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
source('utils/mwas_federated.R');source('utils/mwas_poisson.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
root <- file.path(run_dir,'federated');dir.create(root,recursive=TRUE,showWarnings=FALSE)
reps <- as.integer(Sys.getenv('MWAS_SIM_REPS','200'))
sizes <- as.integer(strsplit(Sys.getenv('MWAS_SIM_SIZES','12,100,600'),',',fixed=TRUE)[[1]])
windows <- as.integer(strsplit(Sys.getenv('MWAS_SIM_WINDOWS','7,28'),',',fixed=TRUE)[[1]])
stopifnot(is.finite(reps),reps>=2,length(sizes)>0,all(is.finite(sizes)&sizes>=10),all(windows %in% c(3,7,14,28)))
input_md5 <- unname(tools::md5sum(file.path(run_dir,'private','matched_exposures.rds')))
protocol_id <- if(file.exists(file.path(root,'protocol.json')))unname(tools::md5sum(file.path(root,'protocol.json'))) else 'unavailable'
seed <- as.integer(Sys.getenv('MWAS_SIM_SEED','20261005'));set.seed(seed)
m <- readRDS(file.path(run_dir,'private','matched_exposures.rds'))
raw <- list();i <- 0L
for(adj in mwas_adjustment_labels)for(ex in c('pm25','o3'))for(w in windows)for(n in sizes) {
 available <- mwas_common_windows(m,ex,weather_adjusted=adj=='weather_adjusted')
 if(uniqueN(available$stratum)<n)next
 ids <- sample(unique(available$stratum),n,replace=FALSE)
 d <- copy(available[stratum %in% ids]);exposure <- paste0(ex,'_lag1_',w)
 design <- tryCatch(mwas_count_design(d,exposure,if(ex=='pm25')5 else 10,weather_adjusted=adj=='weather_adjusted'),error=function(e)NULL)
 if(is.null(design)){message('Unidentifiable simulation design: ',ex,' ',w,' ',n);next}
 sets <- split(seq_len(nrow(d)),d$stratum)
 dates <- seq(min(d$date),max(d$date),by='day');date_index <- match(d$date,dates)
 # Known nuisance slopes on spline coordinates; do not simulate exposure-only
 # truth while claiming adjustment was validated under meteorologic effects.
 nuisance <- rep(c(.15,-.10,.10),length.out=ncol(design$expanded_X)-1L)
 for(scenario in c('independent','shared_calendar_ar1'))for(target_or in c(1,1.5))for(rep in seq_len(reps)) {
  truth <- c(log(target_or),nuisance)
  eta <- drop(design$expanded_X %*% truth)
  if(scenario=='shared_calendar_ar1') {
   shock <- as.numeric(arima.sim(list(ar=.7),n=length(dates),sd=.7))
   eta <- eta+shock[date_index]
  }
  cases <- integer(nrow(d))
  for(ii in sets)cases[ii[sample.int(length(ii),1,prob=exp(eta[ii]-max(eta[ii])))]] <- 1L
  y <- numeric(nrow(design$X));s <- rowsum(cases,design$expanded_map,reorder=TRUE)
  y[as.integer(rownames(s))] <- s[,1]
  f <- mwas_count_fit(design,y,expanded_case=cases)
  for(method in c('patient','quasi','hac')) {
   i <- i+1L;ok <- f$status=='ok';se <- if(ok)unname(f$se[method]) else NA_real_
   beta <- if(ok)f$log_or else NA_real_
   raw[[i]] <- data.table(model_adjustment=adj,pollutant=ex,window_days=w,n_events=n,scenario=scenario,
    target_or=target_or,replicate=rep,inference_method=method,status=f$status,
    log_or=beta,se=se,reject=if(ok)abs(beta/se)>qnorm(.975) else FALSE,
    covered=if(ok)abs(beta-log(target_or))<=qnorm(.975)*se else NA)
  }
  if(rep==reps)message('Calibrated ',adj,' ',ex,' ',w,'d n=',n,' ',scenario,' OR=',target_or)
 }
 # Checkpoint aggregate simulation records: no IDs, ZIPs or dates.
 fwrite(rbindlist(raw),file.path(root,'simulation_replicates.csv'))
}
r <- rbindlist(raw)
summary <- r[,{
 valid <- status=='ok';rate <- mean(reject)
 list(attempted=.N,estimable=sum(valid),fit_failure_rate=mean(!valid),
  rejection_rate=rate,mc_se=sqrt(rate*(1-rate)/.N),
  rejection_ci_low=if(sum(reject)==0)0 else qbeta(.025,sum(reject),.N-sum(reject)+1),
  rejection_ci_high=if(sum(reject)==.N)1 else qbeta(.975,sum(reject)+1,.N-sum(reject)),
  conditional_coverage=if(any(valid))mean(covered[valid]) else NA_real_,
  mean_bias=if(any(valid))mean(log_or[valid]-log(target_or[1])) else NA_real_)
},by=.(model_adjustment,pollutant,window_days,n_events,scenario,target_or,inference_method)]
fwrite(summary,file.path(root,'simulation_calibration.csv'))
write_json(list(seed=seed,replicates_per_cell=reps,sizes=sizes,windows=windows,
 truth='Fixed matched sets and actual local exposure series. Primary truth includes exposure and holiday slopes; weather-sensitivity truth additionally includes spline nuisance slopes. Event dates sampled conditionally. This does not test bias from omitting a true weather confounder.',
 shared_shock='Unmeasured site-wide daily AR1=.7, innovation SD=.7 on log risk; independent of exposure generation',
 interpretation='OR1 rejection is empirical type-I error; OR1.5 rejection is power. Failures count as nonrejections; coverage is conditional on estimability.',
 limitations='Conditional cohort simulation, not a complete population/DAG model. Not BH discovery power. Shared shock marginalizes a conditional effect; bias under that scenario need not be zero.',
 matched_exposures_md5=input_md5,completed_cells=nrow(summary),model_adjustments=mwas_adjustment_labels,expected_cells=length(mwas_adjustment_labels)*length(unique(sizes))*length(unique(windows))*2*2*2*3,
 protocol_id=protocol_id),
 file.path(root,'simulation_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
message('Simulation calibration complete: ',nrow(summary),' cells')
