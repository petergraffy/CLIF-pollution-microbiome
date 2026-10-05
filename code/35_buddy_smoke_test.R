#!/usr/bin/env Rscript
rscript <- file.path(R.home('bin'),'Rscript')
tests <- c('tests/test_site_workflow.R','tests/test_mwas.R','tests/test_mwas_power.R','tests/test_mwas_windows.R',
 'tests/test_mwas_federated.R','tests/test_mwas_poisson.R','tests/test_mwas_characteristics.R','tests/smoke_buddy_pipeline.R')
for(test in tests) {
 status <- system2(rscript,test)
 if(status!=0)stop('Smoke test failed: ',test)
}
message('All R site smoke tests passed; synthetic tests require no patient data or network access.')
