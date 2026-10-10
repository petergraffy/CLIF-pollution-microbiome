# Design-based local normal approximation; does not use fitted organism effects.
# Assumes independent matched sets and equal case probabilities under a joint
# exposure/nuisance null. Simulation is required to validate final eligibility.
mwas_design_information <- function(d, exposure, unit, weather_adjusted=FALSE) {
  d <- data.table::copy(d)
  d[, x := get(exposure)/unit]
  d <- d[is.finite(x) & !is.na(holiday)]
  if(weather_adjusted)d <- d[is.finite(tmean_lag1_7) & is.finite(rhmean_lag1_7)]
  keep <- d[, .(ok=sum(case)==1L && sum(case==0L)>=1L &&
                  diff(range(x))>1e-8), by=stratum][ok==TRUE,stratum]
  d <- d[stratum %in% keep]
  out <- list(n_events=data.table::uniqueN(d$stratum),
              n_patients=data.table::uniqueN(d$patient_id),
              information=NA_real_, status='insufficient_design')
  if (!nrow(d)) return(out)
  z <- tryCatch(if(weather_adjusted) cbind(x=d$x, splines::ns(d$tmean_lag1_7,df=3),
                     splines::ns(d$rhmean_lag1_7,df=3), holiday=d$holiday) else
                     cbind(x=d$x,holiday=d$holiday),error=function(e) NULL)
  if (is.null(z)) return(out)
  # Expected conditional likelihood information at zero coefficients is the
  # sum of within-set covariance matrices, each using uniform weights 1/m.
  groups <- split(seq_len(nrow(d)), d$stratum)
  for (ii in groups) z[ii,] <- sweep(z[ii,,drop=FALSE],2,
                                    colMeans(z[ii,,drop=FALSE]),'-')/sqrt(length(ii))
  residual <- qr.resid(qr(z[,-1,drop=FALSE]),z[,1])
  out$information <- sum(residual^2)
  out$status <- if(out$information>1e-10) 'approximation_available' else 'no_residual_information'
  out
}
mwas_normal_power <- function(log_or, information, alpha) {
  if (!is.finite(information) || information<=0) return(NA_real_)
  delta <- abs(log_or)*sqrt(information)
  critical <- qnorm(1-alpha/2)
  pnorm(delta-critical)+pnorm(-delta-critical)
}
mwas_required_information <- function(log_or, alpha, power=.8) {
  stopifnot(abs(log_or)>0,alpha>0,alpha<1,power>alpha,power<1)
  delta <- uniroot(function(x)pnorm(x-qnorm(1-alpha/2))+
                     pnorm(-x-qnorm(1-alpha/2))-power,c(0,20))$root
  (delta/abs(log_or))^2
}
