#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
flow <- fread(file.path(run_dir,'cohort_flow.csv'))
models <- fread(file.path(run_dir,'primary_mwas_models.csv'))
all_models <- fread(file.path(run_dir,'mwas_models.csv'))
exploratory <- fread(file.path(run_dir,'exploratory_mwas_models_min50.csv'))
clinical <- fread(file.path(run_dir,'clinical_summary.csv'))
severity <- fread(file.path(run_dir,'severity_defined_outcomes_exploratory.csv'))
coverage <- fread(file.path(run_dir,'exposure_coverage.csv'))
ab <- fread(file.path(run_dir,'antibiotic_qc.csv'))
b <- readRDS(file.path(run_dir,'private','cohort.rds'))
e48 <- b$episodes[hours_from_icu<=48]
ab48 <- e48[,.(n_episodes=.N,n_hospitalizations=uniqueN(hospitalization_id)),by=ab_status]
fwrite(ab48,file.path(run_dir,'antibiotic_qc_48h.csv'))
summary <- list(site='UCMC',cohort_flow=flow,clinical_summary=clinical,antibiotic_48h=ab48,
 primary_named_organisms=uniqueN(models$organism),primary_valid_models=sum(models$status=='ok'),
 primary_fdr_hits=models[status=='ok' & q_value<.05],
 severity_defined_fdr_hits=severity[status=='ok' & q_value<.05],
 exploratory_valid_models=sum(exploratory$status=='ok'),exploratory_fdr_hits=exploratory[status=='ok' & q_value<.05],
 primary_top_associations=models[status=='ok'][order(p_value)][1:min(10,sum(models$status=='ok'))],
 primary_family='48h ICU culture window; both pollutants and all named organism hypotheses; minimum 100 admissions',
 secondary_note='50-event threshold, alternative culture windows, antibiotics, inpatient referents, and clinical-group results are exploratory')
write_json(summary,file.path(run_dir,'summary.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
escape <- function(x) {x <- as.character(x);x <- gsub('&','&amp;',x,fixed=TRUE);x <- gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
table_html <- function(d) {
 if(!nrow(d))return('<p>None.</p>')
 paste0('<table><thead><tr>',paste0('<th>',escape(names(d)),'</th>',collapse=''),'</tr></thead><tbody>',
 paste(apply(as.data.frame(d),1,function(r)paste0('<tr>',paste0('<td>',escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table>')
}
ranked <- models[status=='ok'][order(p_value)]
ranked <- ranked[,.(organism,pollutant,n_events,odds_ratio=round(odds_ratio,3),ci_low=round(ci_low,3),ci_high=round(ci_high,3),
 p_value=signif(p_value,3),q_value=signif(q_value,3))]
html <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>UCMC acute respiratory MWAS</title>',
 '<style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#172335}h1,h2{color:#173f53}table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:8px;border-bottom:1px solid #dce3e8;text-align:left}th{background:#edf3f6}img{max-width:100%}.note{background:#fff4d8;padding:16px;border-radius:8px}</style></head><body>',
 '<h1>UCMC acute respiratory organism MWAS</h1><p>First-pass development run • 2018–2024 • PM2.5 and ozone</p>',
 '<p class="note">Exploratory local results. Culture detection does not establish infection. Clinical culture selection, prior treatment, residence-to-ZCTA linkage, and incomplete severity measurements remain limitations. This report is not approved for external release.</p>',
 '<h2>Primary analysis</h2><p>ICU entry within 24 hours of hospital admission; respiratory cultures within 48 hours after first ICU entry. Admission-based time-stratified case-crossover, matched on weekday/month/year and residence. Exposure: days 1–7 before admission. Patient-cluster robust standard errors; temperature, humidity, and federal holiday adjustment. PM2.5 per 5 ug/m3; ozone per 10 ppb.</p>',
 '<p>',uniqueN(models$organism),' eligible organism categories; ',sum(models$status=='ok'),' valid pollutant models; ',sum(models$q_value<.05,na.rm=TRUE),' associations with FDR &lt; 0.05.</p>',
 '<h2>Cohort flow</h2>',table_html(flow),'<h2>Primary estimates</h2>',table_html(ranked),
 '<img src="primary_mwas_forest.png" alt="Primary MWAS forest plot"><img src="primary_mwas_screen.png" alt="Primary organism-wide screen">',
 '<h2>Antibiotics before cultures within 48h</h2>',table_html(ab48),
 '<p>Hospitalization counts overlap across treatment categories when patients had multiple specimens. No prior documented antibacterial treatment does not establish absence of outpatient or outside-hospital therapy. Equal timestamps are uncertain.</p>',
 '<h2>Clinical characterization</h2>',table_html(clinical),
 '<p>Charlson uses present-on-admission diagnoses. Modified SOFA uses first 6h/24h ICU data, measured paired P/F, and creatinine-only renal scoring. Incomplete totals remain missing. Scores are not case-only adjustment covariates in the case-crossover model.</p>',
 '<h2>Exposure linkage diagnostics</h2>',table_html(coverage),
 '<h2>Exploratory analyses</h2><p>Separate files report 24h/72h culture windows, pre-antibacterial specimens, removal of known inpatient referents, a lower 50-admission threshold, and severity-defined presentations. These are not independent validation. Distributed lag/nonlinear models and formal interaction tests remain future work.</p>',
 '<h2>Exploratory clinical-group signals</h2>',table_html(severity[status=='ok' & q_value<.05,.(organism,pollutant,clinical_group,n_events,odds_ratio,ci_low,ci_high,p_value,q_value)]),
 '<p>Subgroup signals require follow-up; significance in one group does not establish effect modification. No formal interaction test has been performed.</p>',
 '<h2>Files</h2><ul><li><a href="primary_mwas_models.csv">Primary models</a></li><li><a href="mwas_models.csv">Main and sensitivity models</a></li>',
 '<li><a href="exploratory_mwas_models_min50.csv">Exploratory minimum-50 screen</a></li><li><a href="severity_defined_outcomes_exploratory.csv">Exploratory clinical-group results</a></li>',
 '<li><a href="manifest.json">Analysis manifest and R session</a></li><li><a href="exposure_manifest.json">Exposure sources and checksums</a></li></ul></body></html>')
writeLines(html,file.path(run_dir,'report.html'))
# Preserve executed source alongside outputs, without copying site configuration.
source_files <- c('code/21_prepare_acute_mwas.R','code/22_cache_mwas_exposures.py','code/23_run_acute_mwas.R','code/24_summarize_acute_mwas.R',
 'utils/mwas.R','utils/mwas_severity.R','resources/mwas/clif_intermittent_med_categories.csv')
dir.create(file.path(run_dir,'source'),showWarnings=FALSE)
file.copy(source_files,file.path(run_dir,'source'),overwrite=TRUE)
write_json(as.list(tools::md5sum(source_files)),file.path(run_dir,'source_checksums.json'),auto_unbox=TRUE,pretty=TRUE)
message('Report: ',file.path(run_dir,'report.html'))
print(ranked)
