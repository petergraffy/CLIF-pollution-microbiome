suppressPackageStartupMessages({library(data.table);library(survival)})
source('utils/mwas.R')
stopifnot(identical(mwas_zip(c('02115','2115','60637-1234','bad','00000')),c('02115','02115','60637',NA_character_,NA_character_)))
r <- mwas_referents(as.Date('2020-02-29'))
stopifnot(length(r)==5,all(lubridate::wday(r)==7),all(format(r,'%Y-%m')=='2020-02'),as.Date('2020-02-29') %in% r)
e <- data.table(episode_id=c('e1','e2'),hospitalization_id='h1',collect_dttm=as.POSIXct(c('2020-01-01 12:00:00','2020-01-01 14:00:00'),tz='UTC'))
m <- data.table(hospitalization_id='h1',admin_dttm=as.POSIXct('2020-01-01 12:00:00',tz='UTC'),med_category='ceftriaxone',med_route_category='iv')
a <- mwas_ab_history(e,m,'h1')
stopifnot(a$ab_status[1]=='timing_uncertain',a$ab_status[2]=='documented_prior',a$hours_since_first_ab[2]==2)
set.seed(72)
d <- data.table(stratum=rep(1:200,each=4),patient_id=rep(1:200,each=4),case=rep(c(1,0,0,0),200),
 pm25=rnorm(800),tmean_lag1_7=rnorm(800),rhmean_lag1_7=rnorm(800),holiday=0L)
f <- mwas_fit(d,'pm25',1,100)
stopifnot(f$n_events==200,f$n_referents==600,f$status=='ok',is.finite(f$p_value))
d[stratum<=25 & case==1,pm25:=NA_real_]
f <- mwas_fit(d,'pm25',1,100)
stopifnot(f$n_events==175,f$n_referents==525)
stopifnot(identical(mwas_negative_text(c('no legionella isolated','coagulase negative staphylococcus species','candida sp. not candida albicans','unable to isolate legionella')),c(TRUE,FALSE,FALSE,TRUE)))
cat('MWAS synthetic tests passed: ZIP normalization, referent calendar, antibiotic ordering, matched-set completeness and model fit.\n')
source('utils/mwas_severity.R')
t0 <- as.POSIXct('2020-01-01 12:00:00',tz='UTC')
co <- data.table(hospitalization_id='h',icu_in=t0)
fixture <- list(
 labs=data.table(hospitalization_id='h',lab_collect_dttm=t0,lab_result_dttm=t0,
  lab_category=c('creatinine','bilirubin_total','platelet_count','po2_arterial'),lab_value_numeric=c(1.5,2,110,80)),
 vitals=data.table(hospitalization_id='h',recorded_dttm=t0,vital_category=c('map','weight'),vital_value=c(65,80)),
 respiratory_support=data.table(hospitalization_id='h',recorded_dttm=t0,device_category='IMV',fio2_set=.5),
 patient_assessments=data.table(hospitalization_id='h',recorded_dttm=t0,assessment_category='gcs_total',numerical_value=13),
 medication_admin_continuous=data.table(hospitalization_id=character(),admin_dttm=as.POSIXct(character()),med_category=character(),
 med_dose=numeric(),med_dose_unit=character(),mar_action_group=character()))
reader <- function(table,cols,...)copy(fixture[[table]])
s <- mwas_severity(co,reader)
stopifnot(s$sofa_24h_total==9,s$sofa_24h_resp==3,s$sofa_24h_observed_domains==6)
fixture$labs <- fixture$labs[lab_category!='bilirubin_total']
s <- mwas_severity(co,reader)
stopifnot(is.na(s$sofa_24h_total),s$sofa_24h_observed_domains==5,s$sofa_24h_partial_sum==7)
# A gas and FiO2 more than one hour apart must not create a respiratory score.
fixture$respiratory_support[,recorded_dttm:=t0+7200]
s <- mwas_severity(co,reader)
stopifnot(is.na(s$sofa_24h_resp))
cat('Severity fixtures passed: known component scores, missingness, and contemporaneous P/F pairing.\n')
