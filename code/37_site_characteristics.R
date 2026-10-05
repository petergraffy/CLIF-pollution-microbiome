#!/usr/bin/env Rscript
# Clinical linkage stays local; exports are admission-level aggregate statistics.
suppressPackageStartupMessages({library(data.table);library(arrow);library(dplyr);library(jsonlite)})
source('utils/clif_io.R');source('utils/mwas.R');source('utils/mwas_federated.R');source('utils/mwas_characteristics.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
root <- file.path(run_dir,'federated');dir.create(root,recursive=TRUE,showWarnings=FALSE)
b <- readRDS(file.path(run_dir,'private','cohort.rds'));manifest <- read_json(file.path(run_dir,'manifest.json'),simplifyVector=TRUE)
read_columns <- function(table,required,optional=character()) {
 ds <- arrow::open_dataset(find_table_path(table));missing <- setdiff(required,names(ds))
 if(length(missing))stop(table,' missing fields: ',paste(missing,collapse=', '))
 as.data.table(ds %>% select(all_of(c(required,intersect(optional,names(ds))))) %>% collect())
}
h <- read_columns('hospitalization',c('hospitalization_id','patient_id','admission_dttm','discharge_dttm','age_at_admission','zipcode_five_digit','discharge_category'),'admission_type_category')
p <- read_columns('patient',c('patient_id','sex_category','race_category'),'ethnicity_category')
a <- read_columns('adt',c('hospitalization_id','in_dttm','out_dttm','location_category'))
for(v in c('admission_dttm','discharge_dttm'))set(h,j=v,value=mwas_ts(h[[v]]))
for(v in c('in_dttm','out_dttm'))set(a,j=v,value=mwas_ts(a[[v]]))
h[,admission_date:=as.Date(admission_dttm,tz=manifest$timezone)]
h <- h[!is.na(admission_date)&admission_date>=as.Date(manifest$start_date)&admission_date<=as.Date(manifest$end_date)]
if(anyDuplicated(h$hospitalization_id))stop('Hospitalizations must be unique')
h <- mwas_attach_descriptors(h,p,a);h[,year:=as.integer(format(admission_date,'%Y'))]
h[,valid_zip:=!is.na(mwas_zip(zipcode_five_digit))]
co_fields <- intersect(c('hospitalization_id','charlson',grep('^sofa_[0-9]+h_(total|observed_domains)$',names(b$cohort),value=TRUE)),names(b$cohort))
h <- merge(h,b$cohort[,..co_fields],by='hospitalization_id',all.x=TRUE,sort=FALSE)
# All descriptions are admissions, including repeat admissions for a patient.
ids48 <- unique(b$episodes[hours_from_icu<=48,hospitalization_id]);pos48 <- unique(b$detections[hours_from_icu<=48,hospitalization_id])
h[,`:=`(early_icu_valid_zip=hospitalization_id %in% b$cohort$hospitalization_id,
 primary_culture48=hospitalization_id %in% ids48,named_organism48=hospitalization_id %in% pos48)]
e48 <- b$episodes[hours_from_icu<=48];det48 <- b$detections[hours_from_icu<=48]
first <- e48[order(collect_dttm,episode_id),.SD[1],by=hospitalization_id][,.(hospitalization_id,first_culture_icu_hours=hours_from_icu,first_culture_antibiotic_status=ab_status)]
h <- merge(h,first,by='hospitalization_id',all.x=TRUE,sort=FALSE)
m <- readRDS(file.path(run_dir,'private','matched_exposures.rds'))
pm <- unique(mwas_common_windows(m,'pm25')$stratum);o3 <- unique(mwas_common_windows(m,'o3')$stratum)
h[,`:=`(matched_pm25=hospitalization_id %in% pm,matched_o3=hospitalization_id %in% o3)]
h[,matched_both:=primary_culture48&matched_pm25&matched_o3]
# Optional IMV capture: absent records remain unknown; observed device records
# with no IMV indicate no documented IMV, not proof of absence of ventilation.
h[,documented_imv:=NA]
rs_path <- find_table_path('respiratory_support',required=FALSE)
if(!is.na(rs_path)) {
 ds <- arrow::open_dataset(rs_path)
 if(all(c('hospitalization_id','recorded_dttm','device_category') %in% names(ds))) {
  hospital_ids <- h$hospitalization_id
  r <- as.data.table(ds %>% select(hospitalization_id,recorded_dttm,device_category) %>% filter(hospitalization_id %in% hospital_ids) %>% collect())
  r[,`:=`(hospitalization_id=as.character(hospitalization_id),recorded_dttm=mwas_ts(recorded_dttm))]
  r <- merge(r,h[,.(hospitalization_id,admission_dttm,discharge_dttm)],by='hospitalization_id')
  r <- r[!is.na(recorded_dttm)&recorded_dttm>=admission_dttm&recorded_dttm<=discharge_dttm&!is.na(device_category)&trimws(device_category)!='']
  observed <- r[,.(documented_imv=any(mwas_clean(device_category)=='imv')),by=hospitalization_id]
  h[observed,on='hospitalization_id',documented_imv:=i.documented_imv]
 }
}
saveRDS(h,file.path(run_dir,'private','site_characteristics.rds'))
cohorts <- list(early_icu_valid_zip=h$early_icu_valid_zip,primary_culture48=h$primary_culture48,
 primary_culture48_matched_both=h$matched_both,named_organism48=h$named_organism48)
num <- c('age_at_admission','charlson','sofa_6h_total','sofa_24h_total','hospital_los_days','icu_los_days','first_culture_icu_hours')
catvars <- c('sex_category','race_category','ethnicity_category','in_hospital_death','hospice_discharge','documented_imv','first_culture_antibiotic_status')
table1 <- rbindlist(lapply(names(cohorts),function(group) {
 d <- h[cohorts[[group]]];out <- mwas_descriptive_rows(d,num,catvars);out[,cohort:=group];out
}))
table1[,site:=clif_site_name];fwrite(table1,file.path(root,'table1_long.csv'))
fmt <- function(x)if(is.finite(x))format(round(x,1),trim=TRUE) else 'NA'
table1[,display:=vapply(seq_len(.N),function(i) {
 r <- .SD[i];if(r$kind=='continuous')paste0(fmt(r$median),' [',fmt(r$q25),', ',fmt(r$q75),']; observed ',r$n_observed,'/',r$n_total) else
 if(r$kind=='categorical')paste0(r$n,' (',if(is.finite(r$pct))format(round(r$pct,2),trim=TRUE) else 'NA','%)') else as.character(r$n)
},character(1))]
labels <- c('Hospital admissions'='Hospital admissions','Unique patients'='Unique patients',
 age_at_admission='Age, years',charlson='POA Charlson index, points',sofa_6h_total='Complete modified SOFA, first 6 ICU hours',
 sofa_24h_total='Complete modified SOFA, first 24 ICU hours',hospital_los_days='Hospital length of stay, days',
 icu_los_days='Total observed ICU length of stay, days',first_culture_icu_hours='First respiratory culture, hours after ICU entry',
 sex_category='Recorded sex',race_category='Recorded race',ethnicity_category='Recorded ethnicity',
 in_hospital_death='In-hospital death',hospice_discharge='Hospice discharge',documented_imv='Documented invasive ventilation',
 first_culture_antibiotic_status='Antibiotic status at first respiratory culture')
 table1[,characteristic:=unname(labels[variable])]
 table1[,category:=fifelse(level=='TRUE','Yes',fifelse(level=='FALSE','No',level))]
 wide <- dcast(table1,site+characteristic+category~cohort,value.var='display',fill='0 (0%)')
 setcolorder(wide,c('site','characteristic','category','primary_culture48','primary_culture48_matched_both','early_icu_valid_zip','named_organism48'))
 fwrite(wide,file.path(root,'table1.csv'))
years <- seq(as.integer(substr(manifest$start_date,1,4)),as.integer(substr(manifest$end_date,1,4)))
annual_groups <- c(list(clif_hospitalization_records=rep(TRUE,nrow(h)),any_icu=h$any_icu,early_icu_all=h$early_icu),cohorts)
annual <- list();category <- list();i <- 0L
for(y in years)for(group in names(annual_groups)) {
 d <- h[year==y & annual_groups[[group]]];N <- nrow(d)
 eps <- e48[hospitalization_id %in% d$hospitalization_id];det <- det48[hospitalization_id %in% d$hospitalization_id]
 known <- !is.na(d$in_hospital_death);rs_known <- !is.na(d$documented_imv)
 row <- data.table(site=clif_site_name,year=y,cohort=group,n_admissions=N,n_patients=uniqueN(d$patient_id),
 n_early_icu=sum(d$early_icu),n_valid_zip=sum(d$valid_zip),n_cultured48=sum(d$primary_culture48),
 culture48_pct=if(N>0)100*sum(d$primary_culture48)/N else NA_real_,n_culture_episodes48=nrow(eps),
 culture_episodes48_per100_admissions=if(N>0)100*nrow(eps)/N else NA_real_,
 n_named_organism48=sum(d$named_organism48),named_detection_pct_among_cultured=if(nrow(eps)>0)100*sum(d$named_organism48)/sum(d$primary_culture48) else NA_real_,
 n_episodes_with_named_organism=uniqueN(det$episode_id),n_distinct_named_taxa=uniqueN(det$organism_category),
 n_death=sum(d$in_hospital_death,na.rm=TRUE),n_mortality_observed=sum(known),n_mortality_missing=sum(!known),
 hospital_mortality_pct=if(any(known))100*mean(d$in_hospital_death[known]) else NA_real_,
 n_hospice=sum(d$hospice_discharge,na.rm=TRUE),n_documented_imv=sum(d$documented_imv,na.rm=TRUE),
 n_resp_support_observed=sum(rs_known),n_resp_support_missing=sum(!rs_known),
 documented_imv_pct=if(any(rs_known))100*mean(d$documented_imv[rs_known]) else NA_real_)
 for(v in num) {
  x <- if(v %in% names(d))as.numeric(d[[v]]) else rep(NA_real_,N);ok <- is.finite(x)
  q <- if(any(ok))quantile(x[ok],c(.25,.5,.75),names=FALSE) else rep(NA_real_,3)
  row[,c(paste0(v,c('_n_observed','_n_missing','_median','_q25','_q75'))):=list(sum(ok),sum(!ok),q[2],q[1],q[3])]
 }
 i <- i+1L;annual[[i]] <- row
 category[[i]] <- mwas_descriptive_rows(d,character(),c('sex_category','race_category','ethnicity_category','admission_type_category','discharge_category','first_culture_antibiotic_status'))[kind=='categorical'][,`:=`(site=clif_site_name,year=y,cohort=group)]
}
fwrite(rbindlist(annual),file.path(root,'site_year_characteristics.csv'))
fwrite(rbindlist(category,fill=TRUE),file.path(root,'site_year_categories.csv'))
eps <- merge(e48,h[,.(hospitalization_id,year)],by='hospitalization_id')
sources <- eps[,.(n_episodes=.N,n_admissions=uniqueN(hospitalization_id)),by=.(year,fluid_category)]
sources[,`:=`(site=clif_site_name,cohort='primary_culture48')]
sources[,pct_episodes:=100*n_episodes/sum(n_episodes),by=year]
fwrite(sources,file.path(root,'site_year_specimen_sources.csv'))
antibiotic_practices <- eps[,.(n_episodes=.N,n_admissions=uniqueN(hospitalization_id)),by=.(year,ab_status)]
antibiotic_practices[,`:=`(site=clif_site_name,cohort='primary_culture48')]
antibiotic_practices[,pct_episodes:=100*n_episodes/sum(n_episodes),by=year]
fwrite(antibiotic_practices,file.path(root,'site_year_antibiotic_practices.csv'))
organisms <- merge(det48,h[,.(hospitalization_id,year)],by='hospitalization_id')[,.(n_admissions=uniqueN(hospitalization_id),n_episodes=uniqueN(episode_id)),by=.(year,organism_category)]
organisms[,site:=clif_site_name];fwrite(organisms,file.path(root,'site_year_organisms.csv'))
write_json(list(site=clif_site_name,unit='hospital admission, including repeat admissions; distinct patients shown separately',
 primary_cohort='early ICU within 24 hospital hours; valid ZIP; any eligible respiratory culture within 48 ICU hours, including cultures without a named organism',
 matched_cohort='primary culture cohort with complete matched sets for BOTH pollutants across all four durations; descriptive only, models remain pollutant-specific',
 year='local hospital admission year; whole hospitalization outcomes assigned to admission year',
 mortality='Expired discharge among completed admissions with known canonical disposition; Hospice reported separately, not counted as death',
 icu_los='union of observed ICU ADT intervals clipped to hospital stay; incomplete or invalid ICU intervals => missing',
 culture='episodes are a timestamp/fluid/method proxy, not a specimen ID; named detection is not infection adjudication',
 charlson_scope='POA Charlson derived for early-ICU admissions with valid ZIP only; not available for all CLIF hospitalization records',
 sofa_scope='complete totals only, derived for admissions with an eligible culture within 72 ICU hours; no assumption of normal unmeasured domains',
 imv='any documented IMV during admission; absent respiratory-support records => unknown, not unventilated',
 denominator='CLIF hospitalization records are the available extract, not necessarily a census of all hospital admissions',
 disclosure='Aggregate exports can contain small cells; review institutional disclosure rules before transfer'),
 file.path(root,'characteristics_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
message('Table 1 and annual site characteristics: ',root)
