suppressPackageStartupMessages(library(data.table))
source('utils/mwas.R');source('utils/mwas_federated.R')
# Known means; every daily lag is required, and day 28 cannot affect 3/7/14.
daily <- as.data.table(as.list(setNames(seq_len(28),paste0('pm25_lag',seq_len(28)))))
averages <- mwas_window_means(daily,'pm25')
stopifnot(averages$pm25_lag1_3==2,averages$pm25_lag1_7==4,
 averages$pm25_lag1_14==7.5,averages$pm25_lag1_28==14.5)
daily[,pm25_lag28:=NA_real_]
averages <- mwas_window_means(daily,'pm25')
stopifnot(is.na(averages$pm25_lag1_28),averages$pm25_lag1_14==7.5)
daily[,pm25_lag2:=Inf]
stopifnot(is.na(mwas_window_means(daily,'pm25')$pm25_lag1_3))
# Missing long-window referent is removed from ALL windows; missing case drops
# its entire set; zero 28-day variation drops a set even if 7-day varies.
d <- data.table(stratum=rep(1:3,each=3),case=rep(c(1L,0L,0L),3),row_id=1:9,
 tmean_lag1_7=10,rhmean_lag1_7=50,holiday=0)
for(w in c(3,7,14,28))d[,paste0('pm25_lag1_',w):=rep(1:3,3)]
d[row_id %in% c(3,4),pm25_lag1_28:=NA_real_]
d[stratum==3,pm25_lag1_28:=1]
shared <- mwas_common_windows(d,'pm25')
stopifnot(identical(shared$row_id,c(1L,2L)),uniqueN(mwas_common_windows(d,'pm25',7)$stratum)==3)
cat('Window checks passed: exact means, complete daily coverage, identical referents, missing-case and constant-exposure exclusions.\n')

# Missing weather must not select the primary analysis population.
d <- data.table(stratum=rep(1:2,each=3),case=rep(c(1L,0L,0L),2),row_id=1:6,
 tmean_lag1_7=c(NA,10,11,12,13,14),rhmean_lag1_7=50,holiday=0)
for(w in c(3,7,14,28))d[,paste0('pm25_lag1_',w):=rep(1:3,2)]
stopifnot(uniqueN(mwas_common_windows(d,'pm25')$stratum)==2,
 uniqueN(mwas_common_windows(d,'pm25',weather_adjusted=TRUE)$stratum)==1)
d[,c('tmean_lag1_7','rhmean_lag1_7'):=NULL]
stopifnot(uniqueN(mwas_common_windows(d,'pm25')$stratum)==2)
cat('Primary matching does not depend on weather availability.\n')
