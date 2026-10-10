#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
source('utils/mwas_federated.R')
source('utils/mwas_power.R')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2)stop('Usage: Rscript code/28_pool_federated_mwas.R OUTPUT_DIR SITE_ESTIMATES.csv [SITE_ESTIMATES.csv ...]')
out_dir <- args[1];dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
d <- rbindlist(lapply(args[-1],fread),use.names=TRUE)
# Old schemas may pool with each other, but never with an explicitly revised specification.
if(!'model_adjustment' %in% names(d))d[,model_adjustment:='legacy_weather_adjusted']
if(anyNA(d$model_adjustment))stop('Mixed legacy and revised adjustment schemas: pooling refused')
if(!'inference_method' %in% names(d))d[,inference_method:='patient_cluster']
required <- c('site','protocol_id','organism','pollutant','analysis','term','family','exposure_unit','log_or','se','status','n_events','n_patients','null_information')
stopifnot(all(required %in% names(d)))
if(uniqueN(d$protocol_id)!=1L)stop('Incompatible protocol/ACS/code versions: pooling refused')
keys <- c('organism','pollutant','analysis','term','family','exposure_unit','culture_window_hours','exposure_window','inference_method','model_adjustment')
if(anyDuplicated(d[,c('site',keys),with=FALSE]))stop('Duplicate site/model rows: pooling refused')
if(any(d$status=='ok' & (!is.finite(d$log_or)|!is.finite(d$se)|d$se<=0)))stop('Invalid successful estimate')
coverage <- d[,.(sites_attempted=uniqueN(site),sites_ok=uniqueN(site[status=='ok'])),by=keys]
fwrite(coverage,file.path(out_dir,'model_coverage.csv'))
is_interaction <- function(x)startsWith(x,'ses:') | startsWith(x,'modifier:')
valid <- d[status=='ok']
if(!nrow(valid))stop('No estimable site models')
pooled <- valid[,c(mwas_meta(log_or,se),list(n_events=sum(n_events),sum_site_patient_counts=sum(n_patients),
 independent_null_information=if(all(is.finite(null_information)))sum(as.numeric(null_information)) else NA_real_)),by=keys]
pooled[,interpretation:=ifelse(k_sites>=2,'multisite_fixed_effect','single_site_only')]
ses_cov <- valid[is_interaction(analysis),{
 a <- .SD[term=='pollution'];b <- .SD[term %in% c('pollution_ses','pollution_modifier')]
 pairs <- merge(a[,.(site,sa=se,cov_main_interaction)],b[,.(site,sb=se)],by='site')
 list(pooled_main_interaction_covariance=sum((1/pairs$sa^2)/sum(1/a$se^2)*
   (1/pairs$sb^2)/sum(1/b$se^2)*pairs$cov_main_interaction))
},by=.(organism,pollutant,analysis,culture_window_hours,exposure_window,inference_method,model_adjustment)]
pooled <- merge(pooled,ses_cov,by=c('organism','pollutant','analysis','culture_window_hours','exposure_window','inference_method','model_adjustment'),all.x=TRUE,sort=FALSE)
# Separate MWAS families; interaction main effects are reference-SES estimates,
# not independent hypothesis families or the primary pooled MWAS.
pooled[,q_value:=NA_real_]
planned <- coverage[!is_interaction(analysis) | term %in% c('pollution_ses','pollution_modifier'),.N,by=.(family,inference_method,model_adjustment)]
pooled[!is_interaction(analysis) | term %in% c('pollution_ses','pollution_modifier'),q_value:={
 target_family <- .BY$family;target_method <- .BY$inference_method;target_adjustment <- .BY$model_adjustment
 n <- planned$N[planned$family==target_family & planned$inference_method==target_method & planned$model_adjustment==target_adjustment]
 p.adjust(p_value,'BH',n=n)
},by=.(family,inference_method,model_adjustment)]
fwrite(pooled,file.path(out_dir,'pooled_mwas.csv'))
primary <- pooled[analysis=='overall' & term=='pollution' & inference_method=='patient_cluster']
power <- rbindlist(lapply(seq_len(nrow(primary)),function(i) {
 r <- primary[i]
 family_size <- planned[family==r$family & inference_method=='patient_cluster' & model_adjustment==r$model_adjustment,N]
 scenarios <- CJ(target_or=c(1.05,1.1,1.2,1.5,2),alpha=c(.05,.05/family_size))
 scenarios[,approximate_power:=vapply(seq_len(.N),function(j)
   mwas_normal_power(log(target_or[j]),r$independent_null_information,alpha[j]),numeric(1))]
 scenarios[,`:=`(organism=r$organism,pollutant=r$pollutant,exposure_window=r$exposure_window,model_adjustment=r$model_adjustment,
   window_days=as.integer(sub('lag1_','',r$exposure_window)),family=r$family,
   k_sites=r$k_sites,
   n_events=r$n_events,independent_null_information=r$independent_null_information)]
 scenarios
}))
fwrite(power,file.path(out_dir,'pooled_power_scenarios.csv'))
# Keep each site's covariance so joint contrasts can be computed later.
fwrite(valid[is_interaction(analysis),.(site,organism,pollutant,analysis,exposure_window,model_adjustment,term,log_or,se,cov_main_interaction)],
 file.path(out_dir,'interaction_site_covariance.csv'))
fwrite(valid[startsWith(analysis,'ses:'),.(site,organism,pollutant,analysis,exposure_window,model_adjustment,term,log_or,se,cov_main_interaction)],file.path(out_dir,'ses_site_covariance.csv'))
# Contrast definitions are fixed in the shared protocol/registry; retain them
# in pooled modifier rows for interpretation without local clinical linkage.
if(all(c('modifier_reference','modifier_comparison') %in% names(d))) {
 defs <- unique(d[startsWith(analysis,'modifier:'),.(analysis,modifier_reference,modifier_comparison)])
 if(anyDuplicated(defs$analysis))stop('Incompatible modifier contrast definitions')
 pooled <- merge(pooled,defs,by='analysis',all.x=TRUE,sort=FALSE)
 fwrite(pooled,file.path(out_dir,'pooled_mwas.csv'))
}
write_json(list(protocol_id=unique(d$protocol_id),sites=sort(unique(d$site)),
 primary='fixed-effect inverse-variance meta-analysis; random DL estimate and heterogeneity descriptive',
 caution='With one site outputs are local estimates; no between-site heterogeneity is inferable. Patient counts may overlap across sites.',
 inputs=lapply(args[-1],function(f)list(file=basename(f),md5=unname(tools::md5sum(f))))),
 file.path(out_dir,'pool_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
message('Pooled ',nrow(pooled),' model terms from ',uniqueN(valid$site),' sites: ',out_dir)
