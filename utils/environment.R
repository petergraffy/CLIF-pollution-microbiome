# No package installation during analysis: restore the pinned environment beforehand.
mwas_check_environment <- function() {
 if(!file.exists('renv.lock'))stop('Run from the repository root')
 if(!requireNamespace('jsonlite',quietly=TRUE))stop('Restore packages first: Rscript -e \'renv::restore(prompt=FALSE)\'')
 lock <- jsonlite::read_json('renv.lock',simplifyVector=FALSE)
 problems <- character()
 for(p in names(lock$Packages)) {
  actual <- if(requireNamespace(p,quietly=TRUE))packageDescription(p)$Version else 'missing'
  if(!identical(actual,lock$Packages[[p]]$Version))problems <- c(problems,paste0(p,': ',actual,'; expected ',lock$Packages[[p]]$Version))
 }
 if(length(problems))stop('Pinned R environment is incomplete. Run renv::restore(prompt=FALSE).\n',paste(problems,collapse='\n'))
 if(as.character(getRversion())!=lock$R$Version)stop('Use R ',lock$R$Version,' for the shared site protocol')
 invisible(TRUE)
}
