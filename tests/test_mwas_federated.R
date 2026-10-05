suppressPackageStartupMessages({library(data.table);library(survival)})
source('utils/mwas_federated.R')
stopifnot(identical(mwas_diagnosis_group(c('J18.9','J44.1','J96.0','I50.9','486','493.9','999'),
 c('ICD10','ICD10','ICD10','ICD10','ICD9','ICD9','bad')),
 c('pneumonia_aspiration','obstructive_airway','other_respiratory','nonrespiratory',
 'pneumonia_aspiration','obstructive_airway','unmapped')))
dx <- data.table(hospitalization_id=c('a','b','c','c','d'),diagnosis_code=c('J18.9','J44.1','J18.9','I50.9','J96.0'),
 diagnosis_code_format='ICD10',diagnosis_primary=1,poa_present=c(1,NA,1,1,0))
g <- mwas_admission_diagnoses(dx,c('a','b','c','d','e'))
stopifnot(identical(g$diagnosis_group,c('pneumonia_aspiration','missing_poa','ambiguous_primary','primary_not_poa','no_primary_diagnosis')))
set.seed(922)
d <- data.table(stratum=rep(1:1600,each=4),patient_id=rep(1:1600,each=4),
 pm=rnorm(6400),ses=rep(rnorm(1600),each=4),tmean_lag1_7=rnorm(6400),
 rhmean_lag1_7=rnorm(6400),holiday=sample(0:1,6400,TRUE),case=0L)
d[,case:=as.integer(seq_len(.N)==sample(seq_len(.N),1,prob=exp(.25*pm+.2*pm*ses))),by=stratum]
f <- mwas_federated_fit(d,'pm',1,TRUE)
stopifnot(all(f$status=='ok'),abs(f$log_or[1]-.25)<.15,abs(f$log_or[2]-.2)<.15,
 all(is.finite(f$cov_main_interaction)),all(f$se>0),f$n_events[1]==1600)
d[,ses:=0]
stopifnot(all(mwas_federated_fit(d,'pm',1,TRUE)$status=='no_ses_variation'))
# Closed-form inverse-variance result: identical site estimates retain beta,
# halve variance, and have zero heterogeneity.
z <- mwas_meta(c(.4,.4),c(.2,.2))
stopifnot(abs(z$log_or-.4)<1e-12,abs(z$se-.2/sqrt(2))<1e-12,z$I2==0,z$tau2==0)
one <- mwas_meta(.4,.2)
stopifnot(one$k_sites==1,is.na(one$I2),one$se==.2)
cat('Federated checks passed: diagnosis missingness, known interaction recovery, covariance, and meta-analysis.\n')
# Exercise the coordinator with aggregate synthetic exports, including rejection
# of duplicate sites and incompatible protocols. No clinical fixtures are used.
tmp <- tempfile('mwas-meta-');dir.create(tmp)
base <- copy(f)
base[,`:=`(site='A',protocol_id='synthetic-v1',organism='synthetic_taxon',pollutant='pm25',
 analysis='ses:poverty_pct',family='ses_primary',exposure_unit='5 ug/m3',culture_window_hours=48L,exposure_window='lag1_7')]
overall <- copy(base[term=='pollution']);overall[,`:=`(analysis='overall',family='primary',null_information=100)]
base <- rbind(base,overall)
duration <- copy(base)
duration[,`:=`(exposure_window='lag1_3',family=ifelse(analysis=='overall','exposure_duration_sensitivity','ses_primary_duration_sensitivity'),
 log_or=log_or+.1,cov_main_interaction=cov_main_interaction*2)]
base <- rbind(base,duration)
file_a <- file.path(tmp,'a.csv');file_b <- file.path(tmp,'b.csv')
fwrite(base,file_a);other <- copy(base);other[,site:='B'];fwrite(other,file_b)
run_pool <- function(files,out)system2(file.path(R.home('bin'),'Rscript'),
 c('code/28_pool_federated_mwas.R',shQuote(out),shQuote(files)),stdout=file.path(tmp,'log'),stderr=file.path(tmp,'log'))
status <- run_pool(c(file_a,file_b),file.path(tmp,'pool'))
if(status!=0)cat(readLines(file.path(tmp,'log')),sep='\n')
stopifnot(status==0)
pooled <- fread(file.path(tmp,'pool','pooled_mwas.csv'))
stopifnot(all(pooled$k_sites==2),all(abs(pooled[analysis=='overall',se]-f$se[1]/sqrt(2))<1e-8),
 abs(pooled[analysis=='ses:poverty_pct' & exposure_window=='lag1_7',pooled_main_interaction_covariance][1]-f$cov_main_interaction[1]/2)<1e-8,
 abs(pooled[analysis=='ses:poverty_pct' & exposure_window=='lag1_3',pooled_main_interaction_covariance][1]-f$cov_main_interaction[1])<1e-8,
 abs(pooled[analysis=='overall' & exposure_window=='lag1_3',log_or]-f$log_or[1]-.1)<1e-8)
stopifnot(run_pool(c(file_a,file_a),file.path(tmp,'duplicate'))!=0)
other[,protocol_id:='incompatible'];fwrite(other,file_b)
stopifnot(run_pool(c(file_a,file_b),file.path(tmp,'mismatch'))!=0)
cat('Coordinator integration checks passed: two-site pooling, covariance, duplicate and protocol rejection.\n')
