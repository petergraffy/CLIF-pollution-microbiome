# Modified SOFA-97: creatinine-only renal domain, measured P/F only, no urine output.
# Incomplete totals stay NA; observed-domain sums are explicitly labeled partial.
mwas_severity <- function(co,read_cols) {
  ids <- co$hospitalization_id
  message('Deriving POA Charlson and early organ dysfunction summaries')
  lab <- read_cols('labs',c('hospitalization_id','lab_collect_dttm','lab_result_dttm','lab_category','lab_value_numeric'),
                    ids,'lab_category',c('creatinine','bilirubin_total','platelet_count','po2_arterial'))
  lab[, `:=`(hospitalization_id=as.character(hospitalization_id),ts=mwas_ts(lab_collect_dttm),value=as.numeric(lab_value_numeric))]
  lab[is.na(ts),ts:=mwas_ts(lab_result_dttm)]
  vit <- read_cols('vitals',c('hospitalization_id','recorded_dttm','vital_category','vital_value'),ids,'vital_category',c('map','weight_kg','weight'))
  vit[, `:=`(hospitalization_id=as.character(hospitalization_id),ts=mwas_ts(recorded_dttm),value=as.numeric(vital_value))]
  rs <- read_cols('respiratory_support',c('hospitalization_id','recorded_dttm','device_category','fio2_set'),ids)
  rs[, `:=`(hospitalization_id=as.character(hospitalization_id),ts=mwas_ts(recorded_dttm),fio2=as.numeric(fio2_set))]
  rs[fio2>1,fio2:=fio2/100];rs <- rs[fio2>=0.21 & fio2<=1]
  gcs <- read_cols('patient_assessments',c('hospitalization_id','recorded_dttm','assessment_category','numerical_value'),ids,'assessment_category',c('gcs_total'))
  gcs[, `:=`(hospitalization_id=as.character(hospitalization_id),ts=mwas_ts(recorded_dttm),value=as.numeric(numerical_value))]
  vaso <- read_cols('medication_admin_continuous',c('hospitalization_id','admin_dttm','med_category','med_dose','med_dose_unit','mar_action_group'),
                   ids,'med_category',c('norepinephrine','epinephrine','dopamine','dobutamine'))
  vaso[, `:=`(hospitalization_id=as.character(hospitalization_id),ts=mwas_ts(admin_dttm),dose=as.numeric(med_dose))]
  tables <- list(lab=lab,vit=vit,rs=rs,gcs=gcs,vaso=vaso)
  for(k in names(tables)) {
   z <- merge(tables[[k]],co[,.(hospitalization_id,icu_in)],by='hospitalization_id')
   z[,dt:=as.numeric(difftime(ts,icu_in,units='hours'))];z <- z[dt>=0 & dt<=24]
   setkey(z,hospitalization_id);tables[[k]] <- z
  }
  score_threshold <- function(x,cut,reverse=FALSE) {
   if(!length(x) || all(!is.finite(x))) return(NA_real_)
   x <- if(reverse) min(x,na.rm=TRUE) else max(x,na.rm=TRUE)
   if(reverse) sum(x<cut) else sum(x>=cut)
  }
  for(w in c(6,24)) {
   scores <- lapply(ids,function(h){
    l <- tables$lab[.(h)][dt<=w];v <- tables$vit[.(h)][dt<=w];r <- tables$rs[.(h)][dt<=w]
    g <- tables$gcs[.(h)][dt<=w];m <- tables$vaso[.(h)][dt<=w]
    renal <- score_threshold(l[lab_category=='creatinine' & value>=0 & value<=30,value],c(1.2,2,3.5,5))
    liver <- score_threshold(l[lab_category=='bilirubin_total' & value>=0 & value<=100,value],c(1.2,2,6,12))
    coag <- score_threshold(l[lab_category=='platelet_count' & value>=0 & value<=3000,value],c(150,100,50,20),TRUE)
    cns <- score_threshold(g[value>=3 & value<=15,value],c(15,13,10,6),TRUE)
    map <- v[vital_category=='map' & value>=20 & value<=250,value]
    cv <- if(length(map)) as.numeric(min(map)<70) else NA_real_
    weight <- v[vital_category %in% c('weight','weight_kg') & value>=20 & value<=300,value]
    wt <- if(length(weight)) median(weight) else NA_real_
    m <- m[mar_action_group=='administered' & dose>0]
    if(nrow(m)) {
     m[,dose_std:=fifelse(med_dose_unit=='mcg/kg/min',dose,fifelse(med_dose_unit=='mcg/min' & is.finite(wt),dose/wt,NA_real_))]
     vals <- m[is.finite(dose_std),fcase(med_category %in% c('norepinephrine','epinephrine') & dose_std>0.1,4,
      med_category %in% c('norepinephrine','epinephrine') & dose_std>0,3,
      med_category=='dopamine' & dose_std>15,4,med_category=='dopamine' & dose_std>5,3,
      med_category=='dopamine' & dose_std>0,2,med_category=='dobutamine' & dose_std>0,2,default=0)]
     if(length(vals)) cv <- max(c(cv,vals),na.rm=TRUE)
     if(any(!is.finite(m$dose_std))) cv <- NA_real_
    }
    # Pair each arterial gas to closest FiO2/support within +/- one hour.
    p <- l[lab_category=='po2_arterial' & value>=20 & value<=800]
    resp_scores <- numeric()
    if(nrow(p) && nrow(r)) for(i in seq_len(nrow(p))) {
     dd <- abs(as.numeric(difftime(r$ts,p$ts[i],units='hours')))
     j <- which.min(dd)
     if(length(j) && dd[j]<=1) {
      pf <- p$value[i]/r$fio2[j];supported <- mwas_clean(r$device_category[j]) %in% c('imv','nippv','cpap')
      resp_scores <- c(resp_scores,if(pf<100 && supported)4 else if(pf<200 && supported)3 else if(pf<300)2 else if(pf<400)1 else 0)
     }
    }
    resp <- if(length(resp_scores)) max(resp_scores) else NA_real_
    s <- c(resp,coag,liver,cv,cns,renal)
    data.table(hospitalization_id=h,resp=resp,coag=coag,liver=liver,cv=cv,cns=cns,renal=renal,
      total=if(all(is.finite(s)))sum(s) else NA_real_,observed_domains=sum(is.finite(s)),partial_sum=sum(s,na.rm=TRUE))
   })
   sc <- rbindlist(scores);setnames(sc,setdiff(names(sc),'hospitalization_id'),paste0('sofa_',w,'h_',setdiff(names(sc),'hospitalization_id')))
   co <- merge(co,sc,by='hospitalization_id',all.x=TRUE)
  }
  co
}
