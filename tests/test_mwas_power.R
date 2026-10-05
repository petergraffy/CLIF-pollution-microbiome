suppressPackageStartupMessages(library(data.table))
source('utils/mwas_power.R')
set.seed(23)
d <- data.table(stratum=rep(1:80,each=4),patient_id=rep(1:80,each=4),
 case=rep(c(1L,0L,0L,0L),80),pm=rnorm(320),tmean_lag1_7=rnorm(320),
 rhmean_lag1_7=rnorm(320),holiday=sample(0:1,320,TRUE))
a <- mwas_design_information(d,'pm',1)
z <- mwas_design_information(d,'pm',.5)
stopifnot(a$n_events==80,abs(z$information/a$information-4)<1e-8)
# Within-set constant exposure contributes no information.
d[,pm:=stratum]
stopifnot(mwas_design_information(d,'pm',1)$n_events==0)
need <- mwas_required_information(log(1.2),.05)
stopifnot(abs(mwas_normal_power(log(1.2),need,.05)-.8)<1e-5,
 abs(mwas_normal_power(0,need,.05)-.05)<1e-10,
 mwas_required_information(log(1.2),.05/214)>need)
cat('Power checks passed: exposure scaling, constant matched exposures, null size, target power, and multiplicity.\n')
