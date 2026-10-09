# No package installation during analysis: restore the pinned environment beforehand.
mwas_check_environment <- function(lockfile='renv.lock') {
 if(!file.exists(lockfile))stop('Run from the repository root')
 if(!requireNamespace('jsonlite',quietly=TRUE))stop('Restore packages first: Rscript -e \'renv::restore(prompt=FALSE)\'')
 lock <- jsonlite::read_json(lockfile,simplifyVector=FALSE)
 problems <- character()
 for(p in names(lock$Packages)) {
  actual <- if(requireNamespace(p,quietly=TRUE))packageDescription(p)$Version else 'missing'
  if(!identical(actual,lock$Packages[[p]]$Version))problems <- c(problems,paste0(p,': ',actual,'; expected ',lock$Packages[[p]]$Version))
 }
 if(length(problems))stop('Pinned R environment is incomplete. Run renv::restore(prompt=FALSE).\n',paste(problems,collapse='\n'))
 runtime <- list(r_version=as.character(getRversion()),reference_r_version=lock$R$Version,
  matches_reference_r=identical(as.character(getRversion()),lock$R$Version),
  platform=R.version$platform,package_versions=as.list(setNames(vapply(names(lock$Packages),
   function(p)packageDescription(p)$Version,character(1)),names(lock$Packages))))
 if(!runtime$matches_reference_r)warning('Running R ',runtime$r_version,
  '; the reference environment used R ',runtime$reference_r_version,
  '. R version differences do not block the pipeline. Run code/35_buddy_smoke_test.R before clinical analysis.',call.=FALSE)
 invisible(runtime)
}
