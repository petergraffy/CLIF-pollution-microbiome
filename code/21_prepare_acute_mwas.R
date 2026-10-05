#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(arrow);library(dplyr);library(data.table);library(jsonlite);library(comorbidity)})
source('utils/clif_io.R'); source('utils/mwas.R')
site_tz <- config_value(config,'site_timezone',env='MWAS_TIMEZONE',default='America/Chicago')
start_date <- as.Date(Sys.getenv('MWAS_START_DATE','2018-01-01'))
end_date <- as.Date(Sys.getenv('MWAS_END_DATE','2024-12-31'))
run_id <- Sys.getenv('MWAS_RUN_ID',paste0(clif_site_name,'_',format(Sys.time(),'%Y%m%d_%H%M%S')))
run_dir <- file.path('output','mwas',run_id); private <- file.path(run_dir,'private')
dir.create(private,recursive=TRUE,showWarnings=FALSE)
read_cols <- function(table,cols,ids=NULL,category=NULL,values=NULL) {
 ds <- arrow::open_dataset(find_table_path(table)) %>% select(all_of(cols))
 if(!is.null(ids)) ds <- ds %>% filter(hospitalization_id %in% ids)
 if(!is.null(category)) ds <- ds %>% filter(.data[[category]] %in% values)
 as.data.table(ds %>% collect())
}
h <- read_cols('hospitalization',c('patient_id','hospitalization_id','admission_dttm','discharge_dttm','age_at_admission','zipcode_five_digit'))
h[, `:=`(patient_id=as.character(patient_id),hospitalization_id=as.character(hospitalization_id),
 admission_dttm=mwas_ts(admission_dttm),discharge_dttm=mwas_ts(discharge_dttm),zip=mwas_zip(zipcode_five_digit))]
stopifnot(!anyDuplicated(h$hospitalization_id))
h[, admission_date:=as.Date(admission_dttm,tz=site_tz)]
all_h <- copy(h)
h <- h[!is.na(admission_date) & admission_date>=start_date & admission_date<=end_date]
a <- read_cols('adt',c('hospitalization_id','in_dttm','out_dttm','location_category'))
a[, `:=`(hospitalization_id=as.character(hospitalization_id),in_dttm=mwas_ts(in_dttm))]
icu <- a[mwas_clean(location_category)=='icu' & !is.na(in_dttm),.(icu_in=min(in_dttm)),by=hospitalization_id]
co <- merge(h,icu,by='hospitalization_id')
co[, hours_to_icu:=as.numeric(difftime(icu_in,admission_dttm,units='hours'))]
co <- co[hours_to_icu>=0 & hours_to_icu<=24]
flow <- data.table(step=c('hospitalizations_in_date_range','with_icu','icu_within_24h','icu_within_24h_valid_zip'),
 n=c(nrow(h),nrow(merge(h,icu,by='hospitalization_id')),nrow(co),nrow(co[!is.na(zip)])))
co <- co[!is.na(zip)]
message('Early ICU cohort: ',nrow(co),' hospitalizations')
c <- read_cols('microbiology_culture',c('hospitalization_id','organism_id','order_dttm','collect_dttm','result_dttm',
 'fluid_category','fluid_name','organism_name','method_category','method_name','organism_category','organism_group','lab_loinc_code'))
c[, hospitalization_id:=as.character(hospitalization_id)]
c <- c[hospitalization_id %in% co$hospitalization_id]
for(v in c('order_dttm','collect_dttm','result_dttm')) set(c,j=v,value=mwas_ts(c[[v]]))
for(v in c('fluid_category','method_category','organism_category','organism_group','method_name')) set(c,j=v,value=mwas_clean(c[[v]]))
fwrite(c[, .N,by=.(fluid_category,method_category)],file.path(run_dir,'specimen_mapping_qc.csv'))
resp <- c('respiratory_tract','respiratory_tract_lower','nasopharynx_upperairway','oropharynx_tongue_oralcavity')
c <- c[method_category=='culture' & fluid_category %in% resp & !is.na(collect_dttm)]
c <- merge(c,co[,.(hospitalization_id,icu_in,discharge_dttm)],by='hospitalization_id')
c[, hours_from_icu:=as.numeric(difftime(collect_dttm,icu_in,units='hours'))]
fwrite(c[,.(n_rows=.N),by=.(timing=fifelse(hours_from_icu<0,'before_icu',fifelse(hours_from_icu<=72,'icu_0_72h','later')))],file.path(run_dir,'culture_timing_qc.csv'))
c <- c[hours_from_icu>=0 & hours_from_icu<=72 & collect_dttm<=discharge_dttm]
c[, episode_id:=mwas_episode_key(.SD)]
# Preserve distinct organisms; remove serial result updates for each organism identity.
setorder(c,episode_id,organism_category,result_dttm)
c <- unique(c,by=c('episode_id','organism_category'),fromLast=TRUE)
episodes <- unique(c[,.(episode_id,hospitalization_id,collect_dttm,fluid_category,hours_from_icu)])
med <- read_cols('medication_admin_intermittent',c('hospitalization_id','admin_dttm','med_category','med_route_category','mar_action_group'))
med[, `:=`(hospitalization_id=as.character(hospitalization_id),med_category=mwas_clean(med_category),
 admin_dttm=mwas_ts(admin_dttm),mar_action_group=mwas_clean(mar_action_group))]
med <- med[hospitalization_id %in% co$hospitalization_id]
observed_hosp <- unique(med$hospitalization_id)
map <- fread('resources/mwas/clif_intermittent_med_categories.csv')
non_bacterial <- c('acyclovir','amphotericin_b','anidulafungin','caspofungin','cidofovir','fluconazole','foscarnet','ganciclovir',
 'isavuconazole','isavuconazonium','itraconazole','maribavir','micafungin','micafungin_posaconazole_voriconazole','nystatin',
 'oseltamivir','peramivir','posaconazole','rezafungin','ribavirin','valacyclovir','valganciclovir','voriconazole')
ab <- setdiff(map[med_group=='CMS_sepsis_qualifying_antibiotics',med_category],non_bacterial)
fwrite(map[med_category %in% ab,.(med_category,description)],file.path(run_dir,'antibacterial_mapping.csv'))
# Record all antibacterials, including gut-directed agents; routes retained for later spectrum sensitivities.
med_ab <- med[mar_action_group=='administered' & med_category %in% ab & !is.na(admin_dttm)]
message('Annotating ',nrow(episodes),' episodes with prior antibiotic administrations')
episodes <- mwas_ab_history(episodes,med_ab,observed_hosp)
c <- merge(c,episodes[,.(episode_id,ab_status)],by='episode_id')
excluded <- c('no_growth','no growth','unknown','other','other_unspecified','mixed_flora','normal_flora','normal_respiratory_flora',
 'mixed_respiratory_flora','respiratory_flora','normal_oral_flora','none','not_identified','unspecified','bacteria_other','fungus_other','gram_negative_rod','gram_positive_cocci','gram_positive_rod','yeast')
fwrite(c[,.(n_episodes=uniqueN(episode_id),n_hospitalizations=uniqueN(hospitalization_id)),by=.(organism_category,organism_group)],file.path(run_dir,'organism_inventory.csv'))
# Explicit named-organism screen. Generic groups/flora retained only in companion outcomes.
c[, negative_text := mwas_negative_text(organism_name)]
fwrite(c[negative_text==TRUE,.(n_episodes=uniqueN(episode_id)),by=.(organism_category,organism_name)],file.path(run_dir,'negative_result_mapping_qc.csv'))
positive <- c[negative_text==FALSE & !organism_group %in% c('no_growth','no growth') & !is.na(organism_category) & organism_category!='' & !organism_category %in% excluded &
 !grepl('no_growth|mixed.*flora|normal.*flora|unspecified|unidentified',organism_category)]
detections <- unique(positive[,.(hospitalization_id,organism_category,organism_group,episode_id,hours_from_icu,ab_status)])
for(w in c(24,48,72)) flow <- rbind(flow,data.table(step=paste0('respiratory_culture_icu_',w,'h'),n=uniqueN(episodes[hours_from_icu<=w,hospitalization_id])))
fwrite(flow,file.path(run_dir,'cohort_flow.csv'))
fwrite(episodes[,.(n_episodes=.N,n_hospitalizations=uniqueN(hospitalization_id)),by=ab_status],file.path(run_dir,'antibiotic_qc.csv'))
# Charlson: POA-only, both ICD versions, missing diagnosis capture remains missing.
dx <- read_cols('hospital_diagnosis',c('hospitalization_id','diagnosis_code','diagnosis_code_format','poa_present'))
dx[,hospitalization_id:=as.character(hospitalization_id)]; dx <- dx[hospitalization_id %in% co$hospitalization_id]
co[, charlson:=NA_real_]; co[hospitalization_id %in% dx$hospitalization_id,charlson:=0]
poa <- dx[poa_present==1 & !is.na(diagnosis_code)]
poa[,code:=gsub('[^A-Za-z0-9]','',toupper(diagnosis_code))]
cc <- list()
for(v in c('9','10')) {
 z <- poa[grepl(v,diagnosis_code_format),.(id=hospitalization_id,code)]
 if(nrow(z)) {
  z <- comorbidity::comorbidity(as.data.frame(z),id='id',code='code',map=paste0('charlson_icd',v,'_quan'),assign0=TRUE)
  cc[[v]] <- data.table(hospitalization_id=z$id,cci=comorbidity::score(z,weights='charlson',assign0=TRUE))
 }
}
if(length(cc)) {cc <- rbindlist(cc)[,.(cci=max(cci)),by=hospitalization_id];co[cc,on='hospitalization_id',charlson:=i.cci]}
saveRDS(list(cohort=co,all_hospitalizations=all_h,episodes=episodes,detections=detections),file.path(private,'cohort_before_severity.rds'))
source('utils/mwas_severity.R')
fwrite(data.table(zip=sort(unique(co$zip))),file.path(private,'required_zips.csv'))
if(Sys.getenv('MWAS_SKIP_SEVERITY','0')=='1') {
 severity <- co[,.(hospitalization_id)]
 for(w in c(6,24)) {
  severity[,paste0('sofa_',w,'h_total'):=NA_real_]
  severity[,paste0('sofa_',w,'h_partial_sum'):=NA_real_]
  severity[,paste0('sofa_',w,'h_observed_domains'):=0L]
 }
 message('Optional SOFA derivation skipped; totals remain missing')
} else if(Sys.getenv('MWAS_REUSE_SEVERITY','0')=='1' && file.exists(file.path(private,'cohort.rds'))) {
 old <- readRDS(file.path(private,'cohort.rds'))$cohort
 anchors <- merge(co[,.(hospitalization_id,icu_in)],old[,.(hospitalization_id,old_icu=icu_in)],by='hospitalization_id')
 stopifnot(nrow(anchors)==nrow(co),all(anchors$icu_in==anchors$old_icu))
 severity <- old
 message('Reusing explicitly requested severity cache with identical ICU anchors')
} else severity <- mwas_severity(co[hospitalization_id %in% episodes$hospitalization_id],read_cols)
severity_cols <- c('hospitalization_id',grep('^sofa_',names(severity),value=TRUE))
co <- merge(co,severity[,..severity_cols],by='hospitalization_id',all.x=TRUE)
saveRDS(list(cohort=co,all_hospitalizations=all_h,episodes=episodes,detections=detections),file.path(private,'cohort.rds'))
# ZIPs are local linkage inputs; never committed or uploaded.
# Required ZIP linkage input written before severity calculation.
write_json(list(run_id=run_id,site=clif_site_name,tables_path=clif_tables_path,timezone=site_tz,
 start_date=as.character(start_date),end_date=as.character(end_date),culture_windows=c(24,48,72),primary_window=48,
 exposure_windows=c(3,7,14,28),primary_exposure_window=7,organism_threshold='all named taxa attempted; no count filter in federated workflow',age_restriction='none',
 sofa_skipped=Sys.getenv('MWAS_SKIP_SEVERITY')=='1',
 antibiotic_scope='documented current hospitalization administrations; outpatient/transfer history unavailable',
 antibiotic_mapping_md5=unname(tools::md5sum('resources/mwas/clif_intermittent_med_categories.csv')),
 episode_definition='hospitalization + order/collection timestamp + fluid + method name + LOINC; no specimen ID',
 git_commit=system('git rev-parse HEAD',intern=TRUE)),file.path(run_dir,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(run_dir,'output/mwas/latest_run.txt')
message('Prepared run: ',run_dir)
