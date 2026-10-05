suppressPackageStartupMessages({library(data.table);library(survival)})
source('utils/mwas.R');source('utils/mwas_federated.R');source('utils/mwas_poisson.R')
set.seed(701)
templates <- rbindlist(lapply(1:12,function(month) {
 dates <- mwas_referents(as.Date(sprintf('2020-%02d-06',month)))
 CJ(zip=c('01001','01002','01003'),date=dates)[,month:=month]
}))
templates[,`:=`(pm=rnorm(.N),tmean_lag1_7=rnorm(.N),rhmean_lag1_7=rnorm(.N),holiday=sample(0:1,.N,TRUE))]
d <- rbindlist(lapply(1:600,function(i) {
 z <- templates[month==sample(1:12,1) & zip==sample(c('01001','01002','01003'),1)]
 z[,`:=`(stratum=as.character(i),patient_id=as.character(ceiling(i/2)),case=0L)]
 z[sample(.N,1,prob=exp(.3*pm)),case:=1L];z
}))
design <- mwas_count_design(d,'pm',1)
count <- mwas_count_fit(design,expanded_case=d$case)
clogit <- mwas_federated_fit(d,'pm',1)
stopifnot(count$status=='ok',clogit$status=='ok',abs(count$log_or-clogit$log_or)<1e-5,
 abs(count$se['patient']-clogit$se)<1e-5,abs(count$log_or-.3)<.15,
 count$se['quasi']>=count$se['model'],count$se['hac']>0)
# Conditional group constants must cancel from both the likelihood and HAC
# efficient scores; this catches failing to project out group intercepts.
shifted <- design
constants <- matrix(rnorm(length(design$N)*ncol(design$X),sd=2),length(design$N))
shifted$X <- design$X+constants[design$g,,drop=FALSE]
shifted$expanded_X <- design$expanded_X+constants[design$g[design$expanded_map],,drop=FALSE]
check_shift <- mwas_count_fit(shifted,expanded_case=d$case)
stopifnot(check_shift$status=='ok',abs(check_shift$log_or-count$log_or)<1e-6,
 max(abs(check_shift$se-count$se))<1e-6)
# Actual calendar gaps must not be compressed into adjacent HAC observations.
score <- matrix(c(1,2,-1,3),2,2)
dates <- as.Date(c('2020-01-01','2020-01-08'))
stopifnot(max(abs(mwas_calendar_meat(score,dates,1)-crossprod(score)))<1e-10)
meat <- mwas_calendar_meat(score,dates,28)
stopifnot(min(eigen(meat,symmetric=TRUE)$values)>-1e-8)
cat('Count likelihood checks passed: identical clogit estimates/patient SEs, known effect, quasi variance and calendar HAC gaps.\n')
