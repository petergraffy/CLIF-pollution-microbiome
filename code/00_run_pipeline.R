#!/usr/bin/env Rscript
# Familiar site entry point; analysis dependencies are restored separately.
if('--help' %in% commandArgs(trailingOnly=TRUE)) {
 cat('Run from the repository root after renv::restore() and site configuration.\n',
     'Fresh run: Rscript code/00_run_pipeline.R\n',
     'Resume: Rscript code/00_run_pipeline.R --resume --run-id EXISTING_RUN_ID\n',
     'Review and return output/runs/<run_id>/; keep output/mwas/ private.\n',sep='')
 quit(save='no',status=0)
}
source('code/34_run_buddy_site.R')
