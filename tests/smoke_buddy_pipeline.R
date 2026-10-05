# Entire local site pipeline on synthetic CLIF/ACS/exposure tables; no network.
suppressPackageStartupMessages({library(data.table);library(arrow);library(jsonlite);library(digest)})
main <- function() {
set.seed(881)
root <- tempfile('clif-buddy-');dir.create(root);tables <- file.path(root,'tables');dir.create(tables)
run_id <- paste0('synthetic_buddy_',Sys.getpid());run_dir <- file.path('output','mwas',run_id)
id <- sprintf('synthetic_h%04d',1:400)
t0 <- as.POSIXct('2020-01-01 12:00:00',tz='UTC')+sample(0:364,400,TRUE)*86400
h <- data.table(patient_id=sprintf('synthetic_p%04d',ceiling((1:400)/2)),hospitalization_id=id,
 admission_dttm=t0,discharge_dttm=t0+7*86400,age_at_admission=sample(25:85,400,TRUE),discharge_category=sample(c('Home','Expired','Hospice','Missing'),400,TRUE),zipcode_five_digit=sample(c('01001','01002','01003'),400,TRUE))
p <- data.table(patient_id=unique(h$patient_id),sex_category=sample(c('Female','Male','Unknown'),200,TRUE),
 race_category=sample(c('White','Black or African American','Asian','Other','Unknown'),200,TRUE),ethnicity_category='Non-Hispanic')
a <- data.table(hospitalization_id=id,in_dttm=t0+3600,out_dttm=t0+7*86400,location_category='icu')
a[1:5,in_dttm:=t0[1:5]+25*3600]
c <- data.table(hospitalization_id=id[1:300],organism_id=paste0('synthetic_o',1:300),
 order_dttm=t0[1:300]+7200,collect_dttm=t0[1:300]+sample(c(3,30,60),300,TRUE)*3600,
 result_dttm=t0[1:300]+4*86400,fluid_category='respiratory_tract',fluid_name='SYNTHETIC RESPIRATORY',
 method_category='culture',method_name='synthetic culture',lab_loinc_code='synthetic',
 organism_category=sample(c('staphylococcus_aureus','candida_albicans','no_growth'),300,TRUE),organism_group='bacteria')
c[,organism_name:=gsub('_',' ',organism_category)]
c[organism_category=='no_growth',organism_group:='no_growth'];c[1:12,fluid_category:='nasopharynx_upperairway']
med <- data.table(hospitalization_id=id,admin_dttm=t0+3600,med_category='acetaminophen',med_route_category='iv',mar_action_group='administered')
med[1:100,med_category:='ceftriaxone']
dx <- data.table(hospitalization_id=id,diagnosis_code=sample(c('J18.9','J44.1','I50.9','J96.0'),400,TRUE),
 diagnosis_code_format='ICD10',diagnosis_primary=1,poa_present=1)
for(n in c('patient','hospitalization','adt','microbiology_culture','medication_admin_intermittent','hospital_diagnosis')) {
 object <- switch(n,patient=p,hospitalization=h,adt=a,microbiology_culture=c,medication_admin_intermittent=med,hospital_diagnosis=dx)
 arrow::write_parquet(object,file.path(tables,paste0('clif_',n,'.parquet')))
}
cache <- file.path(root,'cache');dir.create(cache)
days <- seq(as.Date('2019-01-01'),as.Date('2021-12-31'),by='day')
env <- CJ(zip=c('01001','01002','01003'),date=days)
env[,`:=`(pm25_ug_m3=exp(rnorm(.N,2,.35)),o3_ppb=30+10*sin(as.numeric(date)/365*2*pi)+rnorm(.N,0,3),
 tmean_c=10+15*sin(as.numeric(date)/365*2*pi)+rnorm(.N),rhmean_pct=50+rnorm(.N,0,5),value_source='synthetic',fill_distance_m=0)]
exposure_records <- list()
for(product in c('pm25','o3','weather')) {
 folder <- file.path(cache,product);dir.create(folder)
 cols <- switch(product,pm25=c('zip','date','pm25_ug_m3','value_source','fill_distance_m'),
  o3=c('zip','date','o3_ppb','value_source','fill_distance_m'),weather=c('zip','date','tmean_c','rhmean_pct'))
 for(month in seq(as.Date('2019-12-01'),as.Date('2020-12-01'),by='month')) {
  start <- as.Date(month,origin='1970-01-01');year <- as.integer(format(start,'%Y'));mm <- as.integer(format(start,'%m'))
  relative <- sprintf('%s/%s_%d_%02d.parquet',product,product,year,mm);file <- file.path(cache,relative)
  arrow::write_parquet(env[format(date,'%Y-%m')==format(start,'%Y-%m'),..cols],file)
  exposure_records[[length(exposure_records)+1L]] <- list(product=product,year=year,month=mm,path=relative,
   sha256=digest::digest(file=file,algo='sha256'))
 }
}
write_json(list(scope='national_public_no_clinical_filter',files=exposure_records),file.path(cache,'manifest.json'),auto_unbox=TRUE)
acs_dir <- file.path(root,'acs');dir.create(acs_dir)
fwrite(data.table(zcta=c('01001','01002','01003'),acs_year=2017,poverty_pct=c(10,20,30),
 no_high_school_pct=c(10,15,20),unemployment_pct=c(4,7,10),median_household_income=c(30000,50000,70000),crowding_pct=c(1,2,3)),file.path(acs_dir,'zcta_ses.csv'))
fwrite(data.table(zcta='01001'),file.path(acs_dir,'acs_variables.csv'))
for(table in c('B17001','B15003','B23025','B19013','B25014'))writeLines('{}',file.path(acs_dir,paste0(table,'_metadata.json')))
acs_files <- list.files(acs_dir,full.names=TRUE)
write_json(list(source='synthetic',acs_year=2017,bundled_files_sha256=as.list(setNames(vapply(acs_files,function(f)
 digest::digest(file=f,algo='sha256'),character(1)),basename(acs_files)))),file.path(acs_dir,'manifest.json'),auto_unbox=TRUE)
config_file <- file.path(root,'config.json')
write_json(list(site_name='SYNTHETIC',tables_path=tables,file_type='parquet'),config_file,auto_unbox=TRUE)
Sys.setenv(CLIF_CONFIG_PATH=config_file,MWAS_RUN_ID=run_id,MWAS_RUN_DIR=run_dir,MWAS_EXPOSURE_CACHE=cache,
 MWAS_ACS_DIR=acs_dir,MWAS_SKIP_SEVERITY='1',MWAS_EXPOSURES_ONLY='1',MWAS_START_DATE='2020-01-01',MWAS_END_DATE='2020-12-31',
 MWAS_SIM_REPS='2',MWAS_SIM_SIZES='12,100',MWAS_SIM_WINDOWS='7')
rscript <- file.path(R.home('bin'),'Rscript')
run <- function(script,args=character()) {
 log <- file.path(root,paste0(basename(script),'.log'))
 status <- system2(rscript,c(script,shQuote(args)),stdout=log,stderr=log)
 if(status!=0){cat(readLines(log),sep='\n');stop('Synthetic stage failed: ',script)}
}
old_latest <- if(file.exists('output/mwas/latest_run.txt'))readLines('output/mwas/latest_run.txt') else NULL
on.exit({if(!is.null(old_latest))writeLines(old_latest,'output/mwas/latest_run.txt') else unlink('output/mwas/latest_run.txt');unlink(run_dir,recursive=TRUE);unlink(root,recursive=TRUE)},add=TRUE)
run('code/00_run_pipeline.R')
export_dir <- file.path('output','runs',run_id)
on.exit(unlink(export_dir,recursive=TRUE),add=TRUE)
stopifnot(file.exists(file.path(export_dir,'export_manifest.json')),nrow(fread(file.path(export_dir,'privacy_audit.csv')))==0,
 !any(grepl('private|[.]rds$',list.files(export_dir,recursive=TRUE))))
bundle <- readRDS(file.path(run_dir,'private','cohort.rds'))
stopifnot(nrow(bundle$cohort)==395,all(is.na(bundle$cohort$sofa_24h_total)))
 t1 <- fread(file.path(run_dir,'federated','table1_long.csv'))
 stopifnot(t1[cohort=='primary_culture48' & variable=='Hospital admissions',n]==uniqueN(bundle$episodes[hours_from_icu<=48,hospitalization_id]))
 annual <- fread(file.path(run_dir,'federated','site_year_characteristics.csv'))
 stopifnot(annual[cohort=='clif_hospitalization_records',n_admissions]==400,
 annual[cohort=='early_icu_valid_zip',n_admissions]==395,
 !any(c('patient_id','hospitalization_id','zip','date') %in% names(annual)))
pool <- file.path(run_dir,'federated','pool_local')
 support <- fread(file.path(run_dir,'federated','modifier_case_support.csv'))
 stopifnot(!any(c('patient_id','hospitalization_id','zip','date') %in% names(support)))
 check_n <- merge(support,estimates <- fread(file.path(run_dir,'federated','site_estimates.csv'))[startsWith(analysis,'modifier:') & term=='pollution_modifier',.(organism,pollutant,analysis,exposure_window,n_events)],by=c('organism','pollutant','analysis','exposure_window'))
 stopifnot(all(check_n$n_events_with_modifier==check_n$n_events))
source('utils/mwas_characteristics.R')
est <- fread(file.path(run_dir,'federated','site_estimates.csv'))
stopifnot(!any(c('patient_id','hospitalization_id','zip','date') %in% names(est)),
 all(c('clinical:culture24','clinical:culture72','clinical:no_documented_antibiotics','clinical:pulmonary_sources','clinical:upper_airway_sources','selection_companion') %in% est$analysis),
 all(mwas_modifier_registry()$analysis %in% est$analysis),
 all(c('lag1_3','lag1_7','lag1_14','lag1_28') %in% est$exposure_window))
stopifnot(!any(est[analysis=='modifier:sofa6',status]=='ok'))
 pooled <- fread(file.path(pool,'pooled_mwas.csv'))
 stopifnot(any(startsWith(pooled$analysis,'modifier:')),all(is.na(pooled[startsWith(analysis,'modifier:') & term=='pollution',q_value])))
checks <- fread(file.path(run_dir,'federated','shared_exposure_diagnostics.csv'))
print(head(checks[status=='ok'][order(-abs(point_difference))],6))
stopifnot(any(checks$status=='ok'),max(abs(checks[status=='ok',point_difference]))<1e-4)
if(!is.null(old_latest))writeLines(old_latest,'output/mwas/latest_run.txt')
cat('Synthetic buddy pipeline passed: preflight, six-table cohort, exposures, every analysis family, count checks, pooling, simulations and report.\n')

}
main()
