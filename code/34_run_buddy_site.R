#!/usr/bin/env Rscript
# No external transmission. Resume cached stages with --resume; smoke is separate.
source('utils/config.R')
Sys.setenv(MWAS_EXPOSURE_CACHE=config_value(config,'mwas_exposure_cache',env='MWAS_EXPOSURE_CACHE',default='data/public/exposures'),
 MWAS_ACS_DIR=config_value(config,'mwas_acs_dir',env='MWAS_ACS_DIR',default='data/public/acs/2017'))
args <- commandArgs(trailingOnly=TRUE)
run_id <- Sys.getenv('MWAS_RUN_ID','buddy_first_pass')
run_dir <- file.path('output','mwas',run_id)
resume <- '--resume' %in% args
if(!resume && file.exists(file.path(run_dir,'private','cohort.rds')))stop('Run exists; choose a new MWAS_RUN_ID or use --resume')
Sys.setenv(MWAS_RUN_ID=run_id,MWAS_RUN_DIR=run_dir)
if(!nzchar(Sys.getenv('MWAS_SKIP_SEVERITY')))Sys.setenv(MWAS_SKIP_SEVERITY='1')
rscript <- file.path(R.home('bin'),'Rscript')
run <- function(program,args=character()) {
 code <- system2(program,args=args)
 if(code!=0)stop('Stage failed: ',paste(args,collapse=' '),' (exit ',code,')')
}
run(rscript,'code/33_site_preflight.R')
if(!resume || !file.exists(file.path(run_dir,'private','cohort.rds')))run(rscript,'code/21_prepare_acute_mwas.R')
python <- Sys.getenv('MWAS_PYTHON','python3')
run(python,c('code/22_cache_mwas_exposures.py','--run-dir',shQuote(run_dir),
 '--first-year',substr(Sys.getenv('MWAS_START_DATE','2018-01-01'),1,4),
 '--last-year',substr(Sys.getenv('MWAS_END_DATE','2024-12-31'),1,4)))
run(rscript,'code/23_run_acute_mwas.R')
run(rscript,'code/37_site_characteristics.R')
run(python,'code/26_cache_acs_zcta_ses.py')
run(rscript,'code/27_run_federated_mwas.R')
run(rscript,'code/31_shared_exposure_checks.R')
pool <- file.path(run_dir,'federated','pool_local')
Sys.setenv(MWAS_POOL_DIR=pool)
run(rscript,c('code/28_pool_federated_mwas.R',shQuote(pool),
 shQuote(file.path(run_dir,'federated','site_estimates.csv')),
 shQuote(file.path(run_dir,'federated','shared_exposure_estimates.csv'))))
run(rscript,'code/32_calibrate_mwas.R')
run(rscript,'code/30_audit_duration_candidates.R')
run(rscript,'code/38_audit_modifier_support.R')
run(rscript,'code/29_report_federated_mwas.R')
message('Local buddy run complete. Review aggregate release rules; never share private/ or local caches.')
