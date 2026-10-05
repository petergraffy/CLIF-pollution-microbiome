#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(jsonlite);library(data.table);library(digest)})
source('utils/exports.R')
run_dir <- Sys.getenv('MWAS_RUN_DIR','')
if(!nzchar(run_dir))stop('Set MWAS_RUN_DIR to the completed working run to export')
destination <- Sys.getenv('MWAS_EXPORT_DIR',file.path('output','runs',basename(run_dir)))
mwas_prepare_exports(run_dir,destination)
