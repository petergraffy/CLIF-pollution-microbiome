# Shared admission-level descriptors and predeclared effect modifiers.
suppressPackageStartupMessages(library(data.table))
mwas_demographic_category <- function(x,type) {
 x <- tolower(trimws(as.character(x)));x[is.na(x)|x==''] <- 'missing'
 allowed <- switch(type,sex=c('female','male','unknown'),
  race=c('white','black or african american','asian','american indian or alaska native',
   'native hawaiian or other pacific islander','other','unknown'),
  ethnicity=c('hispanic','non-hispanic','unknown'))
 x[!x %in% c(allowed,'missing')] <- 'unmapped';x
}
mwas_modifier_registry <- function() {
 r <- data.table(indicator=c('age','sex_male','race_black','race_asian','race_aian','race_nhpi','race_other','charlson','sofa6','sofa24'),
  column=c('age_at_admission','sex_category',rep('race_category',5),'charlson','sofa_6h_total','sofa_24h_total'),
  type=c('continuous',rep('contrast',6),rep('continuous',3)),
  center=c(60,rep(NA_real_,6),2,6,6),scale=c(10,rep(NA_real_,6),1,2,2),
  reference=c('age 60 years','female',rep('white',5),'Charlson 2','modified SOFA 6','modified SOFA 6'),
  comparison=c('per 10 years','male','black or african american','asian','american indian or alaska native',
   'native hawaiian or other pacific islander','other','per 1 point','per 2 points','per 2 points'),
  family=c(rep('demographic_modification',7),'clinical_modification','clinical_modification','clinical_modification_24h_sensitivity'))
 r[,analysis:=paste0('modifier:',indicator)];r
}
mwas_modifier_values <- function(d,reg) {
 if(nrow(reg)!=1L)stop('Exactly one modifier registry row is required')
 if(!reg$column %in% names(d))return(rep(NA_real_,nrow(d)))
 v <- d[[reg$column]]
 if(reg$type=='continuous') {
  v <- as.numeric(v);valid <- is.finite(v)
  if(reg$indicator=='age')valid <- valid & v>=0 & v<=120
  if(reg$indicator=='charlson')valid <- valid & v>=0
  if(reg$indicator %in% c('sofa6','sofa24'))valid <- valid & v>=0 & v<=24
  return(ifelse(valid,(v-reg$center)/reg$scale,NA_real_))
 }
 # A separate model compares each listed category with the same reference;
 # unknown/missing and other contrast categories are excluded, never set to 0.
 ifelse(v==reg$reference,0,ifelse(v==reg$comparison,1,NA_real_))
}
mwas_icu_hours <- function(adt,h) {
 a <- copy(adt[tolower(trimws(location_category))=='icu'])
 a <- merge(a,h[,.(hospitalization_id,admission_dttm,discharge_dttm)],by='hospitalization_id')
 if(!nrow(a))return(data.table(hospitalization_id=character(),icu_hours=numeric(),icu_los_complete=logical()))
 a[,.(icu_hours={
   valid <- !is.na(in_dttm)&!is.na(out_dttm)&out_dttm>in_dttm&!is.na(admission_dttm)&!is.na(discharge_dttm)
   if(!all(valid))NA_real_ else {
    start <- pmax(as.numeric(in_dttm),as.numeric(admission_dttm));end <- pmin(as.numeric(out_dttm),as.numeric(discharge_dttm))
    if(any(end<=start))NA_real_ else {
     oo <- order(start,end);start <- start[oo];end <- end[oo]
     total <- 0;left <- start[1];right <- end[1]
     if(length(start)>1)for(i in 2:length(start)) {
      if(start[i]<=right)right <- max(right,end[i]) else {total <- total+right-left;left <- start[i];right <- end[i]}
     }
     (total+right-left)/3600
    }
   }
  }),by=hospitalization_id][,icu_los_complete:=is.finite(icu_hours)]
}
mwas_attach_descriptors <- function(h,p,adt) {
 h <- copy(h);p <- copy(p)
 h[,`:=`(hospitalization_id=as.character(hospitalization_id),patient_id=as.character(patient_id))]
 p[,patient_id:=as.character(patient_id)]
 if(anyDuplicated(p$patient_id))stop('Patient table must have unique patient_id; demographic linkage refused')
 for(v in c('sex_category','race_category','ethnicity_category')) {
  if(!v %in% names(p))p[, (v):=NA_character_]
  p[, (v):=mwas_demographic_category(get(v),sub('_category','',v))]
 }
 h <- merge(h,p[,.(patient_id,sex_category,race_category,ethnicity_category)],by='patient_id',all.x=TRUE,sort=FALSE)
 for(v in c('sex_category','race_category','ethnicity_category'))h[is.na(get(v)),(v):='missing']
 if(!'admission_type_category' %in% names(h))h[,admission_type_category:=NA_character_]
 if(!'discharge_category' %in% names(h))h[,discharge_category:=NA_character_]
 for(v in c('admission_type_category','discharge_category')) {h[, (v):=tolower(trimws(get(v)))];h[is.na(get(v))|get(v)=='',(v):='missing']}
 complete <- !is.na(h$discharge_dttm)&!is.na(h$admission_dttm)&h$discharge_dttm>=h$admission_dttm
 discharge_levels <- c('home','skilled nursing facility (snf)','expired','acute inpatient rehab facility','hospice',
  'long term care hospital (ltach)','acute care hospital','group home','chemical dependency','against medical advice (ama)',
  'assisted living','other','psychiatric hospital','shelter','jail')
 known <- complete & h$discharge_category %in% discharge_levels
 h[,`:=`(hospital_los_days=ifelse(complete,as.numeric(difftime(discharge_dttm,admission_dttm,units='days')),NA_real_),
  in_hospital_death=ifelse(known,discharge_category=='expired',NA),
  hospice_discharge=ifelse(known,discharge_category=='hospice',NA))]
 adt <- copy(adt);adt[,hospitalization_id:=as.character(hospitalization_id)]
 icu <- adt[tolower(trimws(location_category))=='icu' & !is.na(in_dttm),.(first_icu=min(in_dttm)),by=hospitalization_id]
 h <- merge(h,icu,by='hospitalization_id',all.x=TRUE,sort=FALSE)
 h[,`:=`(any_icu=!is.na(first_icu),hours_to_first_icu=as.numeric(difftime(first_icu,admission_dttm,units='hours')))]
 h[,early_icu:=!is.na(hours_to_first_icu)&hours_to_first_icu>=0&hours_to_first_icu<=24]
 h <- merge(h,mwas_icu_hours(adt,h),by='hospitalization_id',all.x=TRUE,sort=FALSE)
 h[any_icu==FALSE,`:=`(icu_hours=0,icu_los_complete=TRUE)]
 h[,icu_los_days:=icu_hours/24];h
}
mwas_descriptive_rows <- function(d,numeric_vars,category_vars) {
 N <- nrow(d);rows <- list(data.table(variable='Hospital admissions',level='',kind='count',n=N,n_total=N,n_observed=N,n_missing=0L,pct=NA_real_,median=NA_real_,q25=NA_real_,q75=NA_real_),
 data.table(variable='Unique patients',level='',kind='count',n=uniqueN(d$patient_id),n_total=N,n_observed=N,n_missing=0L,pct=NA_real_,median=NA_real_,q25=NA_real_,q75=NA_real_))
 for(v in numeric_vars) {
  x <- if(v %in% names(d))as.numeric(d[[v]]) else rep(NA_real_,N);ok <- is.finite(x)
  q <- if(any(ok))quantile(x[ok],c(.25,.5,.75),names=FALSE) else rep(NA_real_,3)
  rows[[length(rows)+1L]] <- data.table(variable=v,level='',kind='continuous',n=sum(ok),n_total=N,n_observed=sum(ok),n_missing=sum(!ok),pct=NA_real_,median=q[2],q25=q[1],q75=q[3])
 }
 for(v in category_vars) {
  x <- if(v %in% names(d))as.character(d[[v]]) else rep(NA_character_,N);x[is.na(x)|x==''] <- 'missing'
  missing <- sum(x %in% c('missing','unknown','unmapped'))
  levels <- sort(unique(x));if(!length(levels))levels <- 'missing'
  for(level in levels)rows[[length(rows)+1L]] <- data.table(variable=v,level=level,kind='categorical',n=sum(x==level),n_total=N,n_observed=N-missing,n_missing=missing,pct=if(N>0)100*sum(x==level)/N else NA_real_,median=NA_real_,q25=NA_real_,q75=NA_real_)
 }
 rbindlist(rows,fill=TRUE)
}
