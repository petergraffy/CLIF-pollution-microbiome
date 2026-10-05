# Offline validation for normal R-only site runs. Rebuilding data is a maintainer task.
mwas_sha256 <- function(path) digest::digest(file=path,algo='sha256')
mwas_validate_public_data <- function(exposure_dir,acs_dir,first_year=2018L,last_year=2024L) {
 check <- function(path,expected) {
  if(!file.exists(path) || dir.exists(path) || !identical(mwas_sha256(path),expected))
   stop('Missing/corrupt bundled public file: ',path,'. Obtain the complete repository data bundle.')
 }
 path <- file.path(exposure_dir,'manifest.json')
 if(!file.exists(path))stop('Missing national exposure manifest: ',path)
 ex <- jsonlite::read_json(path,simplifyVector=FALSE)
 if(!identical(ex$scope,'national_public_no_clinical_filter'))stop('Expected national public inputs, not a site ZIP subset')
 keys <- character()
 for(r in ex$files) {
  if(!r$product %in% c('pm25','o3','weather') || r$month<1 || r$month>12)stop('Invalid exposure file identity')
  expected <- sprintf('%s/%s_%d_%02d.parquet',r$product,r$product,r$year,r$month)
  if(!identical(r$path,expected))stop('Unsafe exposure manifest path')
  keys <- c(keys,sprintf('%s:%d:%d',r$product,r$year,r$month))
  check(file.path(exposure_dir,r$path),r$sha256)
 }
 if(anyDuplicated(keys))stop('Duplicate exposure months')
 needed <- expand.grid(product=c('pm25','o3','weather'),year=seq.int(first_year,last_year),month=1:12)
 required <- c(sprintf('%s:%d:%d',needed$product,needed$year,needed$month),
               sprintf('%s:%d:12',c('pm25','o3','weather'),first_year-1L))
 if(!all(required %in% keys))stop('Bundled exposures do not cover the study period and 28-day lookback')
 ac <- jsonlite::read_json(file.path(acs_dir,'manifest.json'),simplifyVector=FALSE)
 allowed <- c('zcta_ses.csv','acs_variables.csv',paste0(c('B17001','B15003','B23025','B19013','B25014'),'_metadata.json'))
 if(!setequal(names(ac$bundled_files_sha256),allowed))stop('Incomplete or unsafe ACS bundle manifest')
 for(name in names(ac$bundled_files_sha256))check(file.path(acs_dir,name),ac$bundled_files_sha256[[name]])
 list(exposure_manifest_sha256=mwas_sha256(path),acs_indicator_sha256=mwas_sha256(file.path(acs_dir,'zcta_ses.csv')),
      exposure_files=length(ex$files),acs_year=ac$acs_year)
}
