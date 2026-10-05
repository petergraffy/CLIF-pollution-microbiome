suppressPackageStartupMessages({library(data.table);library(survival)})
source('utils/mwas.R');source('utils/mwas_federated.R');source('utils/mwas_characteristics.R')
t <- as.POSIXct('2020-01-01',tz='UTC')
h <- data.table(hospitalization_id=c('a','b','c','d'),patient_id=c('1','1','2','3'),admission_dttm=t,
 discharge_dttm=c(t+86400,t+86400,t+86400,as.POSIXct(NA)),discharge_category=c('Expired','Hospice','Home','Still Admitted'),age_at_admission=c(60,60,80,NA))
p <- data.table(patient_id=c('1','2','3'),sex_category=c('Female','Male','Unknown'),race_category=c('White','Black or African American','noncanonical'))
a <- data.table(hospitalization_id=c('a','a','b'),in_dttm=t+c(3600,7200,3600),out_dttm=t+c(10800,14400,NA),location_category='icu')
z <- mwas_attach_descriptors(h,p,a)
stopifnot(z[hospitalization_id=='a',icu_hours]==3,z[hospitalization_id=='a',in_hospital_death],
 !z[hospitalization_id=='b',in_hospital_death],z[hospitalization_id=='b',hospice_discharge],
 is.na(z[hospitalization_id=='b',icu_hours]),is.na(z[hospitalization_id=='d',in_hospital_death]),
 z[hospitalization_id=='d',race_category]=='unmapped')
stopifnot(inherits(try(mwas_attach_descriptors(h,rbind(p,p[1]),a),silent=TRUE),'try-error'))
r <- mwas_modifier_registry()
stopifnot(identical(mwas_modifier_values(data.table(race_category=c('white','black or african american','unknown','other')),r[indicator=='race_black']),c(0,1,NA_real_,NA_real_)))
stopifnot(identical(mwas_modifier_values(data.table(age_at_admission=c(60,70,NA,125)),r[indicator=='age']),c(0,1,NA_real_,NA_real_)))
stopifnot(is.na(mwas_modifier_values(data.table(sofa_6h_total=NA_real_,sofa_6h_partial_sum=2),r[indicator=='sofa6'])))
desc <- mwas_descriptive_rows(z,c('age_at_admission'),c('race_category','in_hospital_death'))
stopifnot(desc[variable=='Unique patients',n]==3,desc[variable=='age_at_admission',n_missing]==1,
 sum(desc[variable=='race_category',n])==4)
# Known age interaction recovery through exactly the same scaled modifier input.
set.seed(233)
d <- data.table(stratum=rep(1:1800,each=4),patient_id=rep(1:1800,each=4),pm=rnorm(7200),
 age_at_admission=rep(runif(1800,30,90),each=4),tmean_lag1_7=rnorm(7200),rhmean_lag1_7=rnorm(7200),holiday=sample(0:1,7200,TRUE),case=0L)
d[,ses:=mwas_modifier_values(d,r[indicator=='age'])]
d[,case:=as.integer(seq_len(.N)==sample.int(.N,1,prob=exp(.2*pm+.15*pm*ses))),by=stratum]
f <- mwas_federated_fit(d,'pm',1,TRUE)
stopifnot(all(f$status=='ok'),abs(f$log_or[1]-.2)<.12,abs(f$log_or[2]-.15)<.1)
f[term=='pollution_ses',term:='pollution_modifier']
f[,`:=`(site='A',protocol_id='synthetic-modifier',organism='test',pollutant='pm25',analysis='modifier:age',family='demographic_modification',
 exposure_unit='5 ug/m3',culture_window_hours=48L,exposure_window='lag1_7',inference_method='patient_cluster',modifier_reference='age 60 years',modifier_comparison='per 10 years')]
tmp <- tempfile('modifier-pool-');dir.create(tmp);fwrite(f,file.path(tmp,'a.csv'));other <- copy(f);other[,site:='B'];fwrite(other,file.path(tmp,'b.csv'))
status <- system2(file.path(R.home('bin'),'Rscript'),c('code/28_pool_federated_mwas.R',shQuote(file.path(tmp,'pool')),shQuote(file.path(tmp,'a.csv')),shQuote(file.path(tmp,'b.csv'))),stdout=file.path(tmp,'log'),stderr=file.path(tmp,'log'))
if(status!=0)cat(readLines(file.path(tmp,'log')),sep='\n')
stopifnot(status==0)
pooled <- fread(file.path(tmp,'pool','pooled_mwas.csv'))
stopifnot(all(pooled$k_sites==2),is.na(pooled[term=='pollution',q_value]),is.finite(pooled[term=='pollution_modifier',q_value]),
 abs(pooled[term=='pollution_modifier',se]-f[term=='pollution_modifier',se]/sqrt(2))<1e-8,
 abs(pooled[term=='pollution_modifier',pooled_main_interaction_covariance]-f$cov_main_interaction[1]/2)<1e-8)
cat('Characteristics checks passed: overlapping ICU intervals, mortality/hospice, missingness, modifier coding, known interaction recovery and two-site covariance.\n')
