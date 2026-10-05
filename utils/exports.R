# Curated aggregate exports only; never recursively copy a working run.
mwas_prepare_exports <- function(run_dir,destination) {
 federated <- c('site_estimates.csv','shared_exposure_estimates.csv','protocol.json','characteristics_manifest.json',
 'table1.csv','table1_long.csv','site_year_characteristics.csv','site_year_categories.csv','site_year_specimen_sources.csv',
 'site_year_antibiotic_practices.csv','site_year_organisms.csv','modifier_registry.csv','modifier_capture.csv',
 'modifier_case_support.csv','modifier_candidate_diagnostics.csv','duration_candidate_diagnostics.csv',
 'ses_registry.csv','ses_linkage_qc.csv','diagnosis_capture.csv','exposure_window_qc.csv','model_diagnostics.csv',
 'shared_exposure_diagnostics.csv','simulation_calibration.csv','calibration_review.csv',
 'report.html','annual_site_trends.png','duration_sensitivity_forest.png')
 files <- c(file.path('federated',federated),c('cohort_flow.csv','antibiotic_qc.csv','exposure_coverage.csv','exposure_fill_qc.csv'))
 required <- file.path('federated',c('site_estimates.csv','shared_exposure_estimates.csv','protocol.json','table1.csv','report.html'))
 if(!all(file.exists(file.path(run_dir,required))))stop('Run is incomplete: required aggregate results are missing')
 files <- files[file.exists(file.path(run_dir,files))]
 blocked <- '(^|_)(patient_id|hospitalization_id|episode_id|row_id|collect_dttm|admission_dttm|discharge_dttm)$|^(zip|zipcode_five_digit|zcta|date|stratum|tables_path|source_path)$'
 findings <- data.frame(file=character(),field=character(),stringsAsFactors=FALSE)
 walk <- function(x,file) {
  if(is.list(x)) {
   bad <- names(x)[grepl(blocked,names(x),ignore.case=TRUE)]
   if(length(bad))findings <<- rbind(findings,data.frame(file=file,field=bad))
   for(v in x)walk(v,file)
  }
 }
 for(name in files) {
  path <- file.path(run_dir,name)
  if(file.info(path)$isdir || nzchar(Sys.readlink(path)))stop('Unsafe export source: ',name)
  if(endsWith(name,'.csv')) {
   columns <- names(data.table::fread(path,nrows=0))
   bad <- columns[grepl(blocked,columns,ignore.case=TRUE)]
   if(length(bad))findings <- rbind(findings,data.frame(file=name,field=bad))
  }
  if(endsWith(name,'.json'))walk(jsonlite::read_json(path,simplifyVector=FALSE),name)
 }
 if(nrow(findings))stop('Export privacy audit failed: ',paste(paste(findings$file,findings$field,sep=':'),collapse=', '))
 if(dir.exists(destination))stop('Export destination already exists: ',destination,'. Choose a new export folder.')
 parent <- dirname(destination);dir.create(parent,recursive=TRUE,showWarnings=FALSE)
 temporary <- tempfile('export-staging-',tmpdir=parent);dir.create(temporary)
 on.exit(unlink(temporary,recursive=TRUE),add=TRUE)
 for(name in files) {
  target <- file.path(temporary,name);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  if(!file.copy(file.path(run_dir,name),target))stop('Could not copy export: ',name)
 }
 data.table::fwrite(findings,file.path(temporary,'privacy_audit.csv'))
 protocol <- jsonlite::read_json(file.path(run_dir,'federated','protocol.json'))
 estimates <- data.table::fread(file.path(run_dir,'federated','site_estimates.csv'),select='site')
 manifest <- list(run_id=basename(run_dir),site=unique(estimates$site),protocol_version=protocol$version,
  prepared_utc=format(Sys.time(),tz='UTC',usetz=TRUE),files_sha256=as.list(setNames(vapply(files,function(n)
    digest::digest(file=file.path(temporary,n),algo='sha256'),character(1)),files)),
  privacy_audit='No prohibited identifier columns/JSON keys detected in allowlisted structured files; not a comprehensive de-identification guarantee',
  release_status='requires institutional review; small aggregate cells are not suppressed automatically')
 jsonlite::write_json(manifest,file.path(temporary,'export_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
 writeLines(c('Review this aggregate-only folder under your institutional disclosure rules before sharing.',
 'Check small cells, culture/source mapping, coverage, failed fits and calibration diagnostics.',
 'privacy_audit.csv must be empty of findings. The structural audit does not replace human review.',
 'Do not share the working output/mwas run, its private directory, local config or logs.'),file.path(temporary,'README.txt'))
 if(!file.rename(temporary,destination))stop('Could not install completed export folder')
 message('Prepared aggregate exports: ',destination)
 invisible(destination)
}
