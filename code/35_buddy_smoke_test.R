#!/usr/bin/env Rscript
rscript <- file.path(R.home('bin'),'Rscript')
tests <- c('tests/test_mwas.R','tests/test_mwas_power.R','tests/test_mwas_windows.R',
 'tests/test_mwas_federated.R','tests/test_mwas_poisson.R','tests/test_mwas_characteristics.R','tests/smoke_buddy_pipeline.R')
for(test in tests) {
 status <- system2(rscript,test)
 if(status!=0)stop('Smoke test failed: ',test)
}
python <- Sys.getenv('MWAS_PYTHON','python3')
if(system2(python,'tests/test_acs_ses.py')!=0)stop('ACS definition tests failed')
message('All buddy smoke tests passed; synthetic test did not require patient data or network access.')
