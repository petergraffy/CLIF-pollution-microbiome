# Shared MWAS helpers. No clinical files are read when this file is sourced.
suppressPackageStartupMessages({library(data.table); library(lubridate)})
mwas_ts <- function(x) {
  if (inherits(x, 'POSIXt')) return(as.POSIXct(x, tz='UTC'))
  if (is.numeric(x)) return(as.POSIXct(x, origin='1970-01-01', tz='UTC'))
  suppressWarnings(ymd_hms(x, tz='UTC', quiet=TRUE))
}
mwas_zip <- function(x) {
  x <- trimws(as.character(x)); x <- sub('^([0-9]{5})-[0-9]{4}$','\\1',x)
  numeric_zip <- grepl('^[0-9]{1,5}$',x)
  x[numeric_zip] <- sprintf('%05d',as.integer(x[numeric_zip]))
  x[is.na(x) | !grepl('^[0-9]{5}$',x) | x=='00000'] <- NA_character_
  x
}
mwas_clean <- function(x) tolower(trimws(as.character(x)))
mwas_window_means <- function(d, prefixes, windows=c(3L,7L,14L,28L)) {
  d <- copy(d)
  for(prefix in prefixes)for(w in windows) {
    cols <- paste0(prefix,'_lag',seq_len(w))
    stopifnot(all(cols %in% names(d)))
    values <- as.matrix(d[,..cols])
    set(d,j=paste0(prefix,'_lag1_',w),value=ifelse(rowSums(is.finite(values))==w,rowMeans(values),NA_real_))
  }
  d
}
mwas_referents <- function(date) {
  date <- as.Date(date)
  days <- seq(floor_date(date,'month'), ceiling_date(date,'month')-days(1), by='day')
  days[wday(days)==wday(date)]
}
mwas_episode_key <- function(d) {
  # CLIF 2.1 lacks a specimen/test ID. This proxy may collapse simultaneous identical tests.
  paste(d$hospitalization_id, as.numeric(d$collect_dttm), as.numeric(d$order_dttm),
        d$fluid_category,d$method_name,d$lab_loinc_code,sep='|')
}
mwas_ab_history <- function(episodes, meds, observed_hosp) {
  out <- copy(episodes)
  if(nrow(meds)) {
    m <- copy(meds); setkey(m,hospitalization_id,admin_dttm)
    for(i in seq_len(nrow(out))) {
      h <- out$hospitalization_id[i]; t <- out$collect_dttm[i]
      mh <- m[.(h)]
      z <- mh[admin_dttm<t]
      tie <- mh[admin_dttm==t]
      out[i, `:=`(ab_before=nrow(z)>0, ab_equal_time=nrow(tie)>0,
                   ab_admin_count=nrow(z), ab_agents=paste(sort(unique(z$med_category)),collapse=';'),
                   ab_routes=paste(sort(unique(z$med_route_category)),collapse=';'),
                   hours_since_first_ab=if(nrow(z)) as.numeric(difftime(t,min(z$admin_dttm),units='hours')) else NA_real_,
                   hours_since_last_ab=if(nrow(z)) as.numeric(difftime(t,max(z$admin_dttm),units='hours')) else NA_real_)]
    }
  } else out[, `:=`(ab_before=FALSE,ab_equal_time=FALSE,ab_admin_count=0L,ab_agents='',ab_routes='',hours_since_first_ab=NA_real_,hours_since_last_ab=NA_real_)]
  out[, med_record_observed := hospitalization_id %in% observed_hosp]
  out[, ab_status := fifelse(ab_before,'documented_prior',fifelse(ab_equal_time,'timing_uncertain',
                  fifelse(med_record_observed,'none_documented','med_history_unknown')))]
  out
}
mwas_fit <- function(d, exposure, unit, min_events=100L) {
  # Conditional logistic regression: strata are admissions, not organism rows.
  warnings <- character()
  d <- copy(d); d[, x := get(exposure)/unit]
  d <- d[is.finite(x) & is.finite(tmean_lag1_7) & is.finite(rhmean_lag1_7) & !is.na(holiday)]
  keep <- d[, .(ok=sum(case)==1L && sum(case==0L)>=1L && diff(range(x))>1e-8),by=stratum][ok==TRUE,stratum]
  d <- d[stratum %in% keep]
  n <- uniqueN(d$stratum)
  base <- list(n_events=n,n_patients=uniqueN(d$patient_id),n_referents=sum(d$case==0L),
               status='too_few_events',warning='',odds_ratio=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_value=NA_real_)
  if(n<min_events) return(as.data.table(base))
  form <- case ~ x + splines::ns(tmean_lag1_7,df=3) + splines::ns(rhmean_lag1_7,df=3) +
    holiday + strata(stratum) + cluster(patient_id)
  fit <- tryCatch(withCallingHandlers(survival::clogit(form,data=d,method='efron',
    control=survival::coxph.control(iter.max=50)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error')) {base$status <- 'fit_error';base$warning <- conditionMessage(fit);return(as.data.table(base))}
  b <- coef(fit)['x']; se <- sqrt(diag(vcov(fit)))['x']
  base$warning <- paste(unique(warnings),collapse='; ')
  if(length(warnings) || !is.finite(b) || !is.finite(se) || se<=0) {base$status <- 'fit_warning';return(as.data.table(base))}
  base$status <- 'ok';base$odds_ratio <- exp(b);base$ci_low <- exp(b-1.96*se);base$ci_high <- exp(b+1.96*se)
  base$p_value <- 2*pnorm(-abs(b/se))
  as.data.table(base)
}

# Site ETL text overrides are limited to explicit negative/failed results.
# Do not reject 'coagulase negative' or 'not candida albicans': those can be positive organisms.
mwas_negative_text <- function(x) {
 x <- mwas_clean(x)
 !is.na(x) & grepl('^(no |negative|unable to isolate)|overgrown.*unable to isolate',x)
}
