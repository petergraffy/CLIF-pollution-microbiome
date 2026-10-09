#!/usr/bin/env Rscript
# No external transmission. Resume cached stages with --resume; smoke is separate.
source('utils/environment.R');runtime <- mwas_check_environment()
source('utils/config.R')
source('utils/public_data.R')
Sys.setenv(MWAS_EXPOSURE_CACHE=config_value(config,'mwas_exposure_cache',env='MWAS_EXPOSURE_CACHE',default='data/public/exposures'),
 MWAS_ACS_DIR=config_value(config,'mwas_acs_dir',env='MWAS_ACS_DIR',default='data/public/acs/2017'))
args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag) {i <- match(flag,args);if(is.na(i))return('');if(i==length(args))stop('Missing value for ',flag);args[i+1L]}
requested <- get_arg('--run-id')
if(!nzchar(requested))requested <- Sys.getenv('MWAS_RUN_ID','')
if('--resume' %in% args && !nzchar(requested))stop('Resume requires --run-id EXISTING_RUN_ID')
run_id <- if(nzchar(requested))requested else paste0(gsub('[^A-Za-z0-9_-]','_',config$site_name),'_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid())
if(!grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',run_id))stop('Invalid run ID')
run_dir <- file.path('output','mwas',run_id)
resume <- '--resume' %in% args
if(!resume && file.exists(file.path(run_dir,'private','cohort.rds')))stop('Run exists; choose a new MWAS_RUN_ID or use --resume')
Sys.setenv(MWAS_RUN_ID=run_id,MWAS_RUN_DIR=run_dir)
if(!nzchar(Sys.getenv('MWAS_SKIP_SEVERITY')))Sys.setenv(MWAS_SKIP_SEVERITY=if(isTRUE(config$derive_sofa))'0' else '1')
rscript <- file.path(R.home('bin'),'Rscript')
run <- function(program,args=character()) {
 code <- system2(program,args=args)
 if(code!=0)stop('Stage failed: ',paste(args,collapse=' '),' (exit ',code,')')
}
run(rscript,'code/33_site_preflight.R')
public <- mwas_validate_public_data(Sys.getenv('MWAS_EXPOSURE_CACHE'),Sys.getenv('MWAS_ACS_DIR'),
 as.integer(substr(Sys.getenv('MWAS_START_DATE','2018-01-01'),1,4)),as.integer(substr(Sys.getenv('MWAS_END_DATE','2024-12-31'),1,4)))
dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(public,file.path(run_dir,'exposure_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
jsonlite::write_json(runtime,file.path(run_dir,'runtime_environment.json'),pretty=TRUE,auto_unbox=TRUE)
if(!resume || !file.exists(file.path(run_dir,'private','cohort.rds')))run(rscript,'code/21_prepare_acute_mwas.R')
run(rscript,'code/23_run_acute_mwas.R')
run(rscript,'code/37_site_characteristics.R')
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
run(rscript,'code/07_prepare_site_exports.R')
message('Completed run: ',run_id,'; review and return output/runs/',run_id,'/')
