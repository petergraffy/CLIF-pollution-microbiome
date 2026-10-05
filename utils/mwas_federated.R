suppressPackageStartupMessages({library(data.table);library(survival)})

mwas_common_windows <- function(d, exposure, windows=c(3L,7L,14L,28L)) {
  d <- copy(d);cols <- paste0(exposure,'_lag1_',windows)
  if(!all(cols %in% names(d)))stop('Missing exposure windows; rerun code/23_run_acute_mwas.R with MWAS_EXPOSURES_ONLY=1')
  complete <- rowSums(is.finite(as.matrix(d[,..cols])))==length(cols) &
    is.finite(d$tmean_lag1_7) & is.finite(d$rhmean_lag1_7) & !is.na(d$holiday)
  d <- d[complete]
  keep <- d[,.(ok=sum(case)==1L && sum(case==0L)>=1L &&
    all(vapply(.SD,function(x)diff(range(x))>1e-8,logical(1)))),by=stratum,.SDcols=cols][ok==TRUE,stratum]
  d[stratum %in% keep]
}

mwas_diagnosis_group <- function(code, format) {
  code <- gsub('[^A-Z0-9]','',toupper(as.character(code)))
  version <- as.character(format)
  v10 <- grepl('10',version); v9 <- !v10 & grepl('9',version)
  valid <- (v10 & grepl('^[A-Z][0-9]{2}',code)) | (v9 & grepl('^[0-9]{3}',code))
  out <- rep('unmapped',length(code)); out[valid] <- 'nonrespiratory'
  out[valid & ((v10 & grepl('^J',code)) | (v9 & substr(code,1,3)>='460' & substr(code,1,3)<='519'))] <- 'other_respiratory'
  out[valid & ((v10 & grepl('^J(4[1-6])',code)) | (v9 & substr(code,1,3) %in% c('491','492','493','496')))] <- 'obstructive_airway'
  out[valid & ((v10 & grepl('^J(1[2-8]|69)',code)) | (v9 & substr(code,1,3) %in% c(as.character(480:486),'507')))] <- 'pneumonia_aspiration'
  # Influenza with pneumonia is respiratory; retain in other_respiratory unless
  # explicitly included in a later reviewed pneumonia definition.
  out
}
mwas_admission_diagnoses <- function(dx, ids) {
  dx <- copy(dx);dx[,hospitalization_id:=as.character(hospitalization_id)]
  dx <- dx[hospitalization_id %in% ids & diagnosis_primary==1]
  dx[,group:=mwas_diagnosis_group(diagnosis_code,diagnosis_code_format)]
  z <- dx[,.(diagnosis_group=if(any(is.na(poa_present))) 'missing_poa' else
    if(!all(poa_present==1)) 'primary_not_poa' else
    if(uniqueN(group)!=1) 'ambiguous_primary' else group[1]),by=hospitalization_id]
  out <- data.table(hospitalization_id=ids,diagnosis_group='no_primary_diagnosis')
  out[z,on='hospitalization_id',diagnosis_group:=i.diagnosis_group]
  out
}

mwas_federated_fit <- function(d, exposure, unit, interaction=FALSE) {
  d <- copy(d);d[,x:=get(exposure)/unit]
  d <- d[is.finite(x) & is.finite(tmean_lag1_7) & is.finite(rhmean_lag1_7) & !is.na(holiday)]
  if(interaction)d <- d[is.finite(ses)]
  keep <- d[,.(ok=sum(case)==1L && sum(case==0L)>=1L && diff(range(x))>1e-8),by=stratum][ok==TRUE,stratum]
  d <- d[stratum %in% keep]
  terms <- if(interaction)c('pollution','pollution_ses') else 'pollution'
  out <- data.table(term=terms,n_events=uniqueN(d$stratum),n_patients=uniqueN(d$patient_id),
    n_referents=sum(d$case==0),design_rank=NA_integer_,status='insufficient_design',
    log_or=NA_real_,se=NA_real_,cov_main_interaction=NA_real_,null_information=NA_real_)
  if(!nrow(d))return(out)
  if(interaction && uniqueN(d$ses)<2){out[,status:='no_ses_variation'];return(out)}
  # A case-only SES main effect cancels within each stratum; only x:SES is fitted.
  if(interaction)d[,xs:=x*ses]
  form <- if(interaction)case ~ x + xs + splines::ns(tmean_lag1_7,3) + splines::ns(rhmean_lag1_7,3) + holiday + strata(stratum) + cluster(patient_id) else
    case ~ x + splines::ns(tmean_lag1_7,3) + splines::ns(rhmean_lag1_7,3) + holiday + strata(stratum) + cluster(patient_id)
  matrix_form <- if(interaction)~x+xs+splines::ns(tmean_lag1_7,3)+splines::ns(rhmean_lag1_7,3)+holiday else
    ~x+splines::ns(tmean_lag1_7,3)+splines::ns(rhmean_lag1_7,3)+holiday
  z <- tryCatch(model.matrix(matrix_form,d)[,-1,drop=FALSE],error=function(e)NULL)
  if(is.null(z))return(out)
  for(ii in split(seq_len(nrow(d)),d$stratum))z[ii,] <- sweep(z[ii,,drop=FALSE],2,colMeans(z[ii,,drop=FALSE]),'-')/sqrt(length(ii))
  rank <- qr(z)$rank;out[,design_rank:=rank]
  # Degrees-of-freedom safeguard replaces a universal 50/100 event cutoff.
  # It is a computational minimum, not assurance of adequate power or calibration.
  if(out$n_events[1]<=rank+1L || out$n_patients[1]<=rank+1L)return(out)
  core <- if(interaction)c('x','xs') else 'x'
  nuisance <- z[,!colnames(z) %in% core,drop=FALSE]
  residual <- qr.resid(qr(nuisance),z[,core,drop=FALSE])
  info <- crossprod(residual)
  if(min(eigen(info,symmetric=TRUE,only.values=TRUE)$values)<1e-8){out[,status:='unidentifiable_exposure'];return(out)}
  if(!interaction)out[,null_information:=info[1,1]]
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(survival::clogit(form,data=d,method='efron',
    control=survival::coxph.control(iter.max=50)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)NULL)
  if(is.null(fit)){out[,status:='fit_error'];return(out)}
  if(length(warnings)){out[,status:='fit_warning'];return(out)}
  beta <- coef(fit)[core];v <- vcov(fit)[core,core,drop=FALSE]
  if(any(!is.finite(beta)) || any(!is.finite(v)) || min(eigen(v,symmetric=TRUE,only.values=TRUE)$values)<=0){out[,status:='invalid_covariance'];return(out)}
  out[,`:=`(status='ok',log_or=unname(beta),se=sqrt(diag(v)))]
  if(interaction)out[,cov_main_interaction:=v['x','xs']]
  out
}

mwas_meta <- function(beta,se) {
  stopifnot(length(beta)==length(se),length(beta)>0,all(is.finite(beta)),all(is.finite(se)&se>0))
  w <- 1/se^2;k <- length(w);b <- sum(w*beta)/sum(w);s <- sqrt(1/sum(w))
  q <- sum(w*(beta-b)^2);den <- sum(w)-sum(w^2)/sum(w)
  tau2 <- if(k>1 && den>0)max(0,(q-(k-1))/den) else NA_real_
  wr <- if(k>1)1/(se^2+tau2) else w
  br <- sum(wr*beta)/sum(wr);sr <- sqrt(1/sum(wr))
  list(k_sites=k,log_or=b,se=s,odds_ratio=exp(b),ci_low=exp(b-1.96*s),ci_high=exp(b+1.96*s),
    p_value=2*pnorm(-abs(b/s)),Q=if(k>1)q else NA_real_,
    heterogeneity_p=if(k>1)pchisq(q,k-1,lower.tail=FALSE) else NA_real_,
    I2=if(k>1 && q>0)max(0,(q-(k-1))/q)*100 else if(k>1)0 else NA_real_,tau2=tau2,
    random_log_or=br,random_se=sr,random_or=exp(br),random_ci_low=exp(br-1.96*sr),random_ci_high=exp(br+1.96*sr))
}
