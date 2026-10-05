#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(arrow);library(survival);library(ggplot2);library(jsonlite)})
source('utils/mwas.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
cache <- Sys.getenv('MWAS_EXPOSURE_CACHE','data/mwas_cache')
min_events <- as.integer(Sys.getenv('MWAS_MIN_EVENTS','100'))
manifest <- read_json(file.path(run_dir,'manifest.json'),simplifyVector=TRUE)
bundle <- readRDS(file.path(run_dir,'private','cohort.rds'));co <- bundle$cohort
cultured_ids <- unique(bundle$episodes$hospitalization_id)
# Include the entire early-ICU cohort for admission/culture-selection companions.
stopifnot(!anyDuplicated(co$hospitalization_id))
# Federal holiday dates (including observed dates), used identically for case/referent rows.
holidays <- as.Date(character())
for(y in 2017:2025) {
 fixed <- as.Date(paste0(y,c('-01-01','-07-04','-11-11','-12-25')))
 if(y>=2021)fixed <- c(fixed,as.Date(paste0(y,'-06-19')))
 observed <- fixed + ifelse(lubridate::wday(fixed)==7,-1,ifelse(lubridate::wday(fixed)==1,1,0))
 holidays <- c(holidays,fixed,observed)
 for(spec in list(c(1,2,3),c(2,2,3),c(5,2,-1),c(9,2,1),c(10,2,2),c(11,5,4))) {
  ds <- seq(as.Date(sprintf('%d-%02d-01',y,spec[1])),by='day',length.out=31)
  ds <- ds[lubridate::month(ds)==spec[1] & lubridate::wday(ds)==spec[2]]
  holidays <- c(holidays,ds[if(spec[3]<0)length(ds) else spec[3]])
 }
}
matched <- rbindlist(lapply(seq_len(nrow(co)),function(i) {
 h <- co[i];dates <- mwas_referents(h$admission_date)
 data.table(stratum=h$hospitalization_id,patient_id=h$patient_id,zip=h$zip,date=dates,case=as.integer(dates==h$admission_date))
}))
matched[,holiday:=as.integer(date %in% holidays)]
# Flag referents overlapping a known hospitalization; standard matched sets remain primary.
ah <- copy(bundle$all_hospitalizations)
ah[, `:=`(start=as.Date(admission_dttm,tz=manifest$timezone),end=as.Date(discharge_dttm,tz=manifest$timezone))]
ah <- ah[!is.na(start) & !is.na(end),.(patient_id,start,end)]
matched[,row_id:=.I]
blocked <- ah[matched,on=.(patient_id,start<=date,end>=date),nomatch=0,allow.cartesian=TRUE,.(row_id=i.row_id)]
matched[, inpatient_referent:=case==0L & row_id %in% blocked$row_id]
needed_zips <- unique(matched$zip)
read_product <- function(product,cols) {
 files <- list.files(file.path(cache,product),pattern='[.]parquet$',full.names=TRUE)
 if(!length(files))stop('Missing exposure cache: ',product)
 d <- rbindlist(lapply(files,function(f)as.data.table(arrow::read_parquet(f,col_select=tidyselect::all_of(cols)))))
 d[, `:=`(zip=mwas_zip(zip),date=as.Date(date))];d <- d[zip %in% needed_zips]
 stopifnot(!anyDuplicated(d[,.(zip,date)]))
 d
}
message('Loading cached ZCTA exposures')
pm <- read_product('pm25',c('zip','date','pm25_ug_m3','value_source','fill_distance_m'))
o3 <- read_product('o3',c('zip','date','o3_ppb','value_source','fill_distance_m'))
weather <- read_product('weather',c('zip','date','tmean_c','rhmean_pct'))
fill_qc <- rbindlist(list(pm[,.(product='pm25',rows=.N,n_filled=sum(fill_distance_m>0,na.rm=TRUE),max_fill_distance_m=max(fill_distance_m,na.rm=TRUE)),by=.(year=lubridate::year(date))],
 o3[,.(product='o3',rows=.N,n_filled=sum(fill_distance_m>0,na.rm=TRUE),max_fill_distance_m=max(fill_distance_m,na.rm=TRUE)),by=.(year=lubridate::year(date))]))
fwrite(fill_qc,file.path(run_dir,'exposure_fill_qc.csv'))
lookup <- unique(matched[,.(zip,date)]);setkey(lookup,zip,date)
for(prod in c('pm25','o3','tmean','rhmean')) {
 d <- switch(prod,pm25=pm[,.(zip,date,value=pm25_ug_m3)],o3=o3[,.(zip,date,value=o3_ppb)],
  tmean=weather[,.(zip,date,value=tmean_c)],rhmean=weather[,.(zip,date,value=rhmean_pct)])
 setkey(d,zip,date)
 for(lag in 1:28) {
  q <- copy(lookup[,.(zip,date)]);q[, date:=date-lag]
  vals <- d[q,on=.(zip,date),value]
  set(lookup,j=paste0(prod,'_lag',lag),value=vals)
 }
}
lookup <- mwas_window_means(lookup,c('pm25','o3','tmean','rhmean'))
matched <- merge(matched,lookup,by=c('zip','date'),all.x=TRUE)
fwrite(matched[,.(n_rows=.N,n_case_rows=sum(case),n_inpatient_referents=sum(inpatient_referent),
 n_pm25_missing=sum(!is.finite(pm25_lag1_7)),n_o3_missing=sum(!is.finite(o3_lag1_7)),
 n_weather_missing=sum(!is.finite(tmean_lag1_7)|!is.finite(rhmean_lag1_7))),by=.(year=lubridate::year(date))],file.path(run_dir,'exposure_coverage.csv'))
saveRDS(matched,file.path(run_dir,'private','matched_exposures.rds'))
if(Sys.getenv('MWAS_EXPOSURES_ONLY','0')=='1') {
 message('Matched exposure preparation complete; legacy count-filtered models skipped')
 quit(save='no',status=0)
}
# Summaries retain incomplete SOFA totals; partial sums are not labeled full SOFA.
fwrite(co[,.(n_hospitalizations=.N,n_patients=uniqueN(patient_id),median_age=median(age_at_admission,na.rm=TRUE),
 n_under18=sum(age_at_admission<18,na.rm=TRUE),n_charlson_observed=sum(!is.na(charlson)),median_charlson=median(charlson,na.rm=TRUE),
 n_sofa6_complete=sum(!is.na(sofa_6h_total)),n_sofa24_complete=sum(!is.na(sofa_24h_total)),
 median_sofa24_observed_domains=median(sofa_24h_observed_domains,na.rm=TRUE))],file.path(run_dir,'clinical_summary.csv'))
results <- list();idx <- 0L
for(w in c(24,48,72)) {
 e <- bundle$episodes[hours_from_icu<=w]
 det <- bundle$detections[hours_from_icu<=w]
 # Organism group mapping is for plot annotation only.
 orgs <- det[,.(n_hospitalizations=uniqueN(hospitalization_id)),by=organism_category]
 fwrite(orgs,file.path(run_dir,paste0('organism_counts_',w,'h.csv')))
 targets <- c('__any_respiratory_culture','__any_named_organism',orgs[n_hospitalizations>=min_events,organism_category])
 variants <- if(w==48)c('all','pre_antibacterial','no_inpatient_referents') else 'all'
 for(variant in variants)for(org in targets) {
  ids <- if(org=='__any_respiratory_culture')unique(e$hospitalization_id) else if(org=='__any_named_organism')unique(det$hospitalization_id) else unique(det[organism_category==org,hospitalization_id])
  if(variant=='pre_antibacterial') {
   if(org=='__any_respiratory_culture') ids <- unique(e[ab_status=='none_documented',hospitalization_id])
   else if(org=='__any_named_organism')ids <- unique(det[ab_status=='none_documented',hospitalization_id])
   else ids <- unique(det[organism_category==org & ab_status=='none_documented',hospitalization_id])
  }
  d <- matched[stratum %in% ids]
  if(variant=='no_inpatient_referents')d <- d[inpatient_referent==FALSE]
  for(ex in c('pm25','o3')) {
   f <- mwas_fit(d,paste0(ex,'_lag1_7'),if(ex=='pm25')5 else 10,min_events)
   f[, `:=`(organism=org,pollutant=ex,culture_window_hours=w,analysis=variant,
     exposure_window='lag1_7',exposure_unit=if(ex=='pm25')'5 ug/m3' else '10 ppb',n_eligible_events=length(ids))]
   idx <- idx+1L;results[[idx]] <- f
  }
 }
 message('Completed culture window ',w,'h')
}
res <- rbindlist(results,fill=TRUE);res[,q_value:=NA_real_]
# Primary family includes BOTH pollutants and all named organisms at 48h. Companions are separate.
res[!grepl('^__',organism),q_value:=p.adjust(p_value,'BH'),by=.(culture_window_hours,analysis)]
res[grepl('^__',organism),q_value:=p.adjust(p_value,'BH'),by=.(culture_window_hours,analysis)]
fwrite(res,file.path(run_dir,'mwas_models.csv'))
primary <- res[culture_window_hours==48 & analysis=='all' & !grepl('^__',organism)]
fwrite(primary,file.path(run_dir,'primary_mwas_models.csv'))
fwrite(res[,.(models=.N,ok=sum(status=='ok'),warnings=sum(status=='fit_warning'),too_few=sum(status=='too_few_events'),
 fdr_hits=sum(q_value<0.05,na.rm=TRUE)),by=.(culture_window_hours,analysis)],file.path(run_dir,'model_diagnostics.csv'))
p <- primary[status=='ok'];p[,label:=gsub('_',' ',organism)]
if(nrow(p)) {
 p[,label:=factor(label,levels=rev(unique(label[order(p_value)])))]
 g <- ggplot(p,aes(x=odds_ratio,y=label,color=q_value<0.05))+geom_vline(xintercept=1,linetype=2,color='grey60')+
 geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',width=0.2)+geom_point()+
 facet_wrap(~pollutant,scales='free_x')+scale_x_log10()+scale_color_manual(values=c('FALSE'='grey45','TRUE'='#a72b46'),name='FDR < 0.05')+
 labs(x='Odds ratio: PM2.5 per 5 ug/m3; ozone per 10 ppb',y=NULL,title='UCMC acute respiratory organism MWAS',
 subtitle='ICU within 24h; cultures within 48h of ICU entry; exposure days 1–7 before admission')+theme_minimal(base_size=11)
 ggsave(file.path(run_dir,'primary_mwas_forest.png'),g,width=11,height=max(6,uniqueN(p$organism)*0.28),dpi=180)
 g <- ggplot(p,aes(x=reorder(label,-log10(p_value)),y=-log10(p_value),color=q_value<0.05,shape=odds_ratio>1))+
 geom_point(size=2.5)+facet_wrap(~pollutant)+geom_hline(yintercept=-log10(.05),linetype=2)+
 scale_color_manual(values=c('FALSE'='grey45','TRUE'='#a72b46'))+
 labs(x=NULL,y='-log10(p)',color='FDR < 0.05',shape='OR > 1',title='UCMC organism-wide acute pollution screen')+
 theme_minimal()+theme(axis.text.x=element_text(angle=65,hjust=1,size=7))
 ggsave(file.path(run_dir,'primary_mwas_screen.png'),g,width=14,height=7,dpi=180)
}
manifest$model_specification <- 'conditional logistic, Efron (one case per stratum), patient-cluster robust SE; ns(7-day mean temperature,3)+ns(7-day mean RH,3)+federal holiday'
manifest$min_events <- min_events;manifest$fdr_family <- 'all named organisms and both pollutants within each culture-window/analysis; 48h all is primary'
manifest$primary_inpatient_referents <- 'standard weekday-month-year referents retained; separate sensitivity excludes known inpatient dates'
manifest$sofa_definition <- 'modified SOFA-97, creatinine-only renal, measured paired P/F within 1h; incomplete total NA; 6h and 24h windows post ICU entry'
manifest$session_info <- capture.output(sessionInfo())
write_json(manifest,file.path(run_dir,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
message('MWAS completed: ',run_dir)
# Lower-count exploratory screen: separately corrected; does not alter the primary family.
exploratory <- list();idx <- 0L
det48 <- bundle$detections[hours_from_icu<=48]
org50 <- det48[,.(n=uniqueN(hospitalization_id)),by=organism_category][n>=50,organism_category]
for(org in org50)for(ex in c('pm25','o3')) {
 ids <- unique(det48[organism_category==org,hospitalization_id])
 f <- mwas_fit(matched[stratum %in% ids],paste0(ex,'_lag1_7'),if(ex=='pm25')5 else 10,50L)
 f[, `:=`(organism=org,pollutant=ex,culture_window_hours=48L,analysis='exploratory_min50',n_eligible_events=length(ids))]
 idx <- idx+1L;exploratory[[idx]] <- f
}
if(length(exploratory)) {
 er <- rbindlist(exploratory);er[,q_value:=p.adjust(p_value,'BH')]
 fwrite(er,file.path(run_dir,'exploratory_mwas_models_min50.csv'))
}
# Exploratory severity-defined outcomes; SOFA is post-admission and not a baseline confounder.
severity_results <- list();idx <- 0L
for(group in c('charlson_0_2','charlson_3plus','sofa24_0_5','sofa24_6plus')) {
 selected <- switch(group,charlson_0_2=co[!is.na(charlson) & charlson<=2,hospitalization_id],
  charlson_3plus=co[!is.na(charlson) & charlson>=3,hospitalization_id],
  sofa24_0_5=co[!is.na(sofa_24h_total) & sofa_24h_total<=5,hospitalization_id],
  sofa24_6plus=co[!is.na(sofa_24h_total) & sofa_24h_total>=6,hospitalization_id])
 for(org in org50)for(ex in c('pm25','o3')) {
  ids <- intersect(unique(det48[organism_category==org,hospitalization_id]),selected)
  f <- mwas_fit(matched[stratum %in% ids],paste0(ex,'_lag1_7'),if(ex=='pm25')5 else 10,50L)
  f[, `:=`(organism=org,pollutant=ex,clinical_group=group,n_eligible_events=length(ids))]
  idx <- idx+1L;severity_results[[idx]] <- f
 }
}
if(length(severity_results)) {
 sr <- rbindlist(severity_results);sr[,q_value:=p.adjust(p_value,'BH')]
 fwrite(sr,file.path(run_dir,'severity_defined_outcomes_exploratory.csv'))
}
message('Exploratory lower-count and clinical-group screens completed')
