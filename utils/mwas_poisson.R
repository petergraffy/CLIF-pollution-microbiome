# Conditional count likelihood plus separate variance checks. Calendar HAC sums
# scores over all ZCTAs on each date; Bartlett weights use actual calendar lags.
suppressPackageStartupMessages(library(data.table))
mwas_count_design <- function(d,exposure,unit,weather_adjusted=FALSE) {
 d <- copy(d);d[,x:=get(exposure)/unit]
 stopifnot(all(is.finite(d$x)),!anyNA(d$holiday))
 if(weather_adjusted)stopifnot(all(is.finite(d$tmean_lag1_7)),all(is.finite(d$rhmean_lag1_7)))
 X <- model.matrix(mwas_model_formula(weather_adjusted=weather_adjusted,conditional=FALSE),d)[,-1,drop=FALSE]
 # Preserve exact candidate-date sets, including exclusions from missing data.
 signatures <- d[,.(signature=paste(sort(as.character(date)),collapse=',')),by=stratum]
 d[signatures,on='stratum',signature:=i.signature]
 d[,group_key:=paste(zip,signature,sep='|')]
 centered <- copy(X)
 for(ii in split(seq_len(nrow(d)),d$stratum))centered[ii,] <- sweep(X[ii,,drop=FALSE],2,colMeans(X[ii,,drop=FALSE]),'-')
 q <- qr(centered,tol=1e-9);columns <- sort(q$pivot[seq_len(q$rank)])
 if(!1L %in% columns)stop('Exposure is not identifiable')
 X <- X[,columns,drop=FALSE]
 z <- cbind(d[,.(group_key,date,case,patient_id,stratum)],as.data.table(X))
 xn <- colnames(X)
 # Date/ZCTA rows in a collapsed stratum must have identical design values.
 check <- z[,lapply(.SD,function(v)diff(range(v))),by=.(group_key,date),.SDcols=xn]
 if(any(as.matrix(check[,..xn])>1e-8))stop('Cannot aggregate inconsistent shared exposures')
 a <- z[,c(list(y=sum(case)),lapply(.SD,function(v)v[1])),by=.(group_key,date),.SDcols=xn]
 setorder(a,group_key,date)
 g <- match(a$group_key,unique(a$group_key));totals <- as.numeric(rowsum(a$y,g,reorder=FALSE))
 if(any(totals<=0))stop('Conditional groups require positive totals')
 row_map <- match(paste(z$group_key,z$date),paste(a$group_key,a$date))
 list(X=as.matrix(a[,..xn]),g=g,N=totals,y=a$y,date=as.Date(a$date),
  expanded_X=X,expanded_map=row_map,patient=d$patient_id,stratum=d$stratum,
  columns=xn,n_events=sum(d$case),n_patients=uniqueN(d$patient_id),n_referents=sum(d$case==0))
}
mwas_count_components <- function(design,beta,y=design$y) {
 X <- design$X;g <- design$g;N <- design$N
 eta <- drop(X %*% beta);maximum <- as.numeric(tapply(eta,g,max))
 e <- exp(eta-maximum[g]);den <- as.numeric(rowsum(e,g,reorder=FALSE));p <- e/den[g]
 mu <- N[g]*p
 means <- rowsum(X*p,g,reorder=FALSE)
 information <- crossprod(X,X*mu)-crossprod(means,means*N)
 list(nll=sum(N*(log(den)+maximum))-sum(y*eta),gradient=drop(crossprod(X,mu-y)),
  information=information,p=p,mu=mu,means=means)
}
mwas_calendar_meat <- function(score,dates,lag=28L) {
 days <- seq(min(dates),max(dates),by='day')
 S <- matrix(0,length(days),ncol(score))
 index <- match(as.Date(dates),days)
 summed <- rowsum(score,index,reorder=TRUE)
 S[as.integer(rownames(summed)),] <- summed
 meat <- crossprod(S)
 for(k in seq_len(min(lag,nrow(S)-1L))) {
  cross <- crossprod(S[(k+1L):nrow(S),,drop=FALSE],S[seq_len(nrow(S)-k),,drop=FALSE])
  meat <- meat+(1-k/(lag+1))*(cross+t(cross))
 }
 meat
}
mwas_count_fit <- function(design,y=design$y,hac_lag=28L,expanded_case=NULL) {
 p <- ncol(design$X)
 result <- list(status='insufficient_design')
 if(sum(y)<=p+1 || length(y)-length(design$N)-p<=0)return(result)
 opt <- tryCatch(optim(rep(0,p),function(b)mwas_count_components(design,b,y)$nll,
  function(b)mwas_count_components(design,b,y)$gradient,method='BFGS',
  control=list(maxit=150,reltol=1e-10)),error=function(e)NULL)
 if(is.null(opt)||opt$convergence!=0){result$status <- 'fit_error';return(result)}
 c <- mwas_count_components(design,opt$par,y)
 # Refine a converged BFGS fit for agreement with the equivalent clogit MLE.
 for(step in seq_len(8L)) {
  if(max(abs(c$gradient))<1e-8 || rcond(c$information)<1e-10)break
  delta <- tryCatch(solve(c$information,c$gradient),error=function(e)NULL)
  if(is.null(delta))break
  candidate <- opt$par-delta
  next_c <- mwas_count_components(design,candidate,y)
  if(!is.finite(next_c$nll) || next_c$nll>c$nll+1e-8)break
  opt$par <- candidate;c <- next_c
 }
 if(max(abs(c$gradient))>1e-4 || rcond(c$information)<1e-10){result$status <- 'unstable_information';return(result)}
 inverse <- tryCatch(solve(c$information),error=function(e)NULL)
 if(is.null(inverse)){result$status <- 'invalid_covariance';return(result)}
 residual_df <- length(y)-length(design$N)-p
 phi <- sum((y-c$mu)^2/pmax(c$mu,1e-12))/residual_df
 # Efficient scores account for the conditioned group intercepts before
 # splitting groups across calendar clusters. Uncentered scores are incorrect.
 score <- (design$X-c$means[design$g,,drop=FALSE])*(y-c$mu)
 hac <- inverse %*% mwas_calendar_meat(score,design$date,hac_lag) %*% inverse
 variance <- c(model=inverse[1,1],quasi=max(1,phi)*inverse[1,1],hac=hac[1,1])
 if(!is.null(expanded_case)) {
  expanded_score <- design$expanded_X*(expanded_case-c$p[design$expanded_map])
  patient_score <- rowsum(expanded_score,design$patient,reorder=FALSE)
  patient_var <- inverse %*% crossprod(patient_score) %*% inverse
  variance <- c(variance,patient=patient_var[1,1])
 }
 if(any(!is.finite(variance)|variance<=0)){result$status <- 'invalid_covariance';return(result)}
 list(status='ok',log_or=opt$par[1],se=sqrt(variance),dispersion=phi,
  n_conditional_groups=length(design$N),n_calendar_days=length(unique(design$date)),
  n_parameters=p,residual_df=residual_df)
}
