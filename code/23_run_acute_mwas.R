#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(arrow);library(jsonlite);library(dplyr)})
source('utils/mwas.R');source('utils/config.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
cache <- config_value(config,'mwas_exposure_cache',env='MWAS_EXPOSURE_CACHE',default='data/public/exposures')
manifest <- read_json(file.path(run_dir,'manifest.json'),simplifyVector=TRUE)
bundle <- readRDS(file.path(run_dir,'private','cohort.rds'));co <- bundle$cohort
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
 d <- as.data.table(arrow::open_dataset(files) %>% dplyr::select(tidyselect::all_of(cols)) %>%
   dplyr::filter(zip %in% needed_zips) %>% dplyr::collect())
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
message('Matched exposure preparation complete; run code/27_run_federated_mwas.R for organism models')
