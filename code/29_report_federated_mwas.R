#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(ggplot2)})
run_dir <- Sys.getenv('MWAS_RUN_DIR',readLines('output/mwas/latest_run.txt',warn=FALSE)[1])
root <- file.path(run_dir,'federated')
pool <- Sys.getenv('MWAS_POOL_DIR',file.path(root,'pool_ucmc'))
d <- fread(file.path(root,'site_estimates.csv'));p <- fread(file.path(pool,'pooled_mwas.csv'))
for(v in intersect(c('odds_ratio','ci_low','ci_high','q_value','log_or','se'),names(p)))set(p,j=v,value=as.numeric(p[[v]]))
if(!'inference_method' %in% names(p))p[,inference_method:='patient_cluster']
primary <- p[inference_method=='patient_cluster' & model_adjustment=='primary_no_weather']
weather_sensitivity <- p[inference_method=='patient_cluster' & model_adjustment=='weather_adjusted' & analysis=='overall']
esc <- function(x) {x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
table <- function(x)paste0('<table><thead><tr>',paste0('<th>',esc(names(x)),'</th>',collapse=''),
 '</tr></thead><tbody>',paste(vapply(seq_len(nrow(x)),function(i)paste0('<tr>',paste0('<td>',
 esc(unlist(x[i],use.names=FALSE)),'</td>',collapse=''),'</tr>'),character(1)),collapse=''),'</tbody></table>')
summary <- p[!is.na(q_value),.(estimated_terms=.N,FDR_hits=sum(q_value<.05)),by=.(family,inference_method,model_adjustment)]
freq <- d[analysis=='overall' & status=='ok',.(n=max(n_events)),by=organism][order(-n)]
common <- head(freq$organism,6)
effects <- primary[exposure_window=='lag1_7' & organism %in% common & (analysis=='overall' | (analysis=='ses:poverty_pct' & term=='pollution_ses')),
 .(organism,pollutant,analysis,n_events,OR=round(odds_ratio,3),CI_low=round(ci_low,3),CI_high=round(ci_high,3),q=round(q_value,3))]
setorder(effects,analysis,organism,pollutant)
duration <- primary[analysis=='overall' & organism %in% common,
 .(organism,pollutant,days=as.integer(sub('lag1_','',exposure_window)),n_events,
 OR=round(odds_ratio,3),CI_low=round(ci_low,3),CI_high=round(ci_high,3),q=round(q_value,3))]
setorder(duration,organism,pollutant,days)
plot_data <- primary[analysis=='overall' & organism %in% common]
short_labels <- setNames(vapply(strsplit(gsub('_',' ',common),' '),function(x)
 paste0(toupper(substr(x[1],1,1)),'. ',paste(x[-1],collapse=' ')),character(1)),gsub('_',' ',common))
plot_data[,`:=`(days=factor(as.integer(sub('lag1_','',exposure_window)),levels=c(3,7,14,28)),
 organism=factor(gsub('_',' ',organism),levels=gsub('_',' ',common)),
 pollutant=factor(pollutant,levels=c('pm25','o3'),labels=c('PM2.5 per 5 ug/m3','Ozone per 10 ppb')))]
if(nrow(plot_data)) {
 g <- ggplot(plot_data,aes(days,odds_ratio,ymin=ci_low,ymax=ci_high))+
  geom_hline(yintercept=1,linetype='dashed',color='grey60')+
  geom_pointrange(aes(color=exposure_window=='lag1_7'))+
  facet_grid(organism~pollutant,scales='free_y',labeller=labeller(organism=short_labels))+scale_y_log10()+
  scale_color_manual(values=c('FALSE'='#666666','TRUE'='#1769aa'),guide='none')+
  labs(x='Days before admission (seven-day primary in blue)',y='Odds ratio and 95% confidence interval',
   title='Exposure duration sensitivity: six most frequent estimable organisms')+
  theme_bw(base_size=11)
 ggsave(file.path(root,'duration_sensitivity_forest.png'),g,width=10,height=10,dpi=180)
}
site_label <- esc(paste(sort(unique(d$site)),collapse=', '))
audit_path <- file.path(root,'duration_candidate_diagnostics.csv')
audit_html <- if(file.exists(audit_path))paste0('<h2>Candidate stability diagnostics</h2>',
 '<p class="note">An FDR-passing result is a screening candidate, not a validated finding. Compare patient-cluster and model-based uncertainty and leave-one-patient-out failures. A very large odds ratio supported by few admissions can reflect sparse-model instability.</p>',
 table(fread(audit_path))) else ''
optional_table <- function(name)if(file.exists(file.path(root,name)))table(fread(file.path(root,name))) else '<p>Unavailable in this run.</p>'
shared_html <- '<p>Unavailable in this run.</p>'
shared_path <- file.path(root,'shared_exposure_diagnostics.csv')
if(file.exists(shared_path)) {
 checks <- fread(shared_path)
 shared_html <- paste0(table(checks[,.(models=.N),by=.(analysis,model_adjustment,status)]),
  table(checks[analysis=='overall' & exposure_window=='lag1_7' & organism %in% common & status=='ok']))
}
calibration_html <- optional_table('simulation_calibration.csv')
calibration_path <- file.path(root,'simulation_calibration.csv')
if(file.exists(calibration_path)) {
 sim <- fread(calibration_path)
 review <- sim[target_or==1,.(null_cells=.N,
  minimum_null_rejection=min(rejection_rate),maximum_null_rejection=max(rejection_rate),
  cells_with_evidence_of_inflation=sum(rejection_ci_low>.05),
  maximum_fit_failure=max(fit_failure_rate)),by=.(inference_method,model_adjustment)]
 fwrite(review,file.path(root,'calibration_review.csv'))
 flag <- if(any(review$cells_with_evidence_of_inflation>0))'<p class="note">Calibration found null scenarios with rejection rates exceeding 5% beyond Monte Carlo uncertainty. These methods are not uniformly calibrated; technical fit success and FDR correction do not establish valid scientific inference. Review the affected scenarios before confirmatory use.</p>' else '<p>No null scenario showed clear excess rejection beyond Monte Carlo uncertainty; this does not establish universal calibration.</p>'
 sim_manifest_path <- file.path(root,'simulation_manifest.json')
 if(file.exists(sim_manifest_path) && file.exists(file.path(root,'protocol.json'))) {
  sm <- read_json(sim_manifest_path,simplifyVector=TRUE)
  if(!identical(sm$protocol_id,unname(tools::md5sum(file.path(root,'protocol.json')))))
   flag <- paste0(flag,'<p>The displayed baseline calibration comes from an earlier protocol run and is retained as a development diagnostic. It does not validate the new modifier interactions.</p>')
 }
 calibration_html <- paste0(flag,table(review),calibration_html)
}
modifier_effects <- primary[startsWith(analysis,'modifier:') & term=='pollution_modifier',
 .(organism,pollutant,analysis,modifier_reference,modifier_comparison,n_events,
 ratio_of_pollution_ORs=round(odds_ratio,3),CI_low=round(ci_low,3),CI_high=round(ci_high,3),q=round(q_value,3))]
annual_html <- '<p>Unavailable in this run.</p>'
annual_path <- file.path(root,'site_year_characteristics.csv')
if(file.exists(annual_path)) {
 yearly <- fread(annual_path)
 trend <- yearly[cohort=='early_icu_valid_zip',.(year,measure='Cultured within 48 ICU hours (% early ICU admissions)',value=culture48_pct)]
 trend <- rbind(trend,yearly[cohort=='primary_culture48',.(year,measure='Named organism detected (% cultured admissions)',value=named_detection_pct_among_cultured)],
  yearly[cohort=='primary_culture48',.(year,measure='In-hospital death (% known primary-cohort outcomes)',value=hospital_mortality_pct)])
 if(any(is.finite(trend$value))) {
  graph <- ggplot(trend,aes(year,value))+geom_line(color='#1769aa')+geom_point(color='#1769aa')+
   facet_wrap(~measure,ncol=1,scales='free_y')+scale_x_continuous(breaks=sort(unique(trend$year)))+
   labs(x='Hospital admission year',y='Percent',title='Annual site culture practices and outcomes',
    caption='Descriptive rates with different denominators; culture detection is not adjudicated infection.')+theme_bw(base_size=11)
  ggsave(file.path(root,'annual_site_trends.png'),graph,width=10,height=8,dpi=180)
 }
 annual_html <- paste0('<img src="annual_site_trends.png" alt="Annual culture practices and outcomes" style="width:100%;height:auto">',table(yearly[cohort %in% c('early_icu_valid_zip','primary_culture48'),
 .(year,cohort,n_admissions,n_patients,n_cultured48,culture48_pct,n_culture_episodes48,n_named_organism48,
 named_detection_pct_among_cultured,n_death,n_mortality_observed,hospital_mortality_pct,hospital_los_days_median,icu_los_days_median)]))
}
clinical <- primary[startsWith(analysis,'clinical:') | analysis=='selection_companion',.(analysis,organism,pollutant,exposure_window,n_events,OR=odds_ratio,CI_low=ci_low,CI_high=ci_high,q=q_value)]
html <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>Federated respiratory MWAS</title>',
 '<style>body{font:16px system-ui;max-width:1500px;margin:40px auto;padding:0 20px;color:#203040}table{display:block;overflow-x:auto;border-collapse:collapse;width:100%;font-size:14px;margin:20px 0}th,td{padding:8px;text-align:left;border-bottom:1px solid #ddd}th{background:#edf2f6}h2{margin-top:32px}.note{background:#fff4d9;padding:16px}</style></head><body>',
 '<h1>Federated respiratory culture MWAS: ',site_label,' development run</h1>',
 '<p class="note">This run contains one site. Results reproduce local estimates and test the aggregation workflow; they are not multisite evidence. Passing model diagnostics does not establish adequate power or reliable sparse-sample inference.</p>',
 '<p>Hospital admissions entering ICU within 24 hours; cultures in the first 48 ICU hours. Mean PM2.5 and ozone exposures on days 1–3, 1–7 (primary), 1–14, and 1–28 before hospital admission. Four CLIF respiratory source categories are included. Models match weekday within month/year. Primary models adjust for holidays without weather terms. Separate weather sensitivities add lag-1–7 temperature and relative-humidity splines (3 df each) across all windows and analysis families. Primary eligibility does not require weather; sensitivity eligibility additionally requires complete weather.</p>',
 '<h2>Table 1: primary cohort and context cohorts</h2>',optional_table('table1.csv'),
 '<p>The primary cohort includes any eligible respiratory culture within 48 ICU hours, including cultures without a named organism. Summaries count hospital admissions; unique patients are shown separately. Median [IQR] rows show observed/total counts. The matched-both cohort is descriptive; pollutant models retain their own matched sets. Six-/24-hour modified SOFA totals require all six domains and are derived only for admissions cultured within 72 hours. Partial scores are not substituted.</p>',
 '<h2>Annual site characteristics</h2>',annual_html,
 '<p>Year is local hospital admission year. CLIF admission counts describe available records, not necessarily the hospital admission census. In-hospital death is Expired discharge; hospice is separate. Mortality percentages use completed admissions with known disposition. ICU duration merges overlapping ADT intervals; incomplete intervals are missing. Full annual exports include category distributions, specimen sources, antibiotic practices and named-organism counts. These are descriptive changes in case mix and documentation, not pollution-effect estimates.</p>',
 '<h2>Clinical and demographic effect modification</h2>',optional_table('modifier_capture.csv'),table(modifier_effects),
 '<h3>Interaction candidate stability</h3>',optional_table('modifier_candidate_diagnostics.csv'),
 '<p>Model-level category/score support is exported in modifier_case_support.csv. A large total event count can hide a contrast supported by only one or two events in a category. FDR-passing candidates remain exploratory when support, model-based uncertainty or leave-one-patient-out refits indicate instability.</p>',
 '<p>Seven-day models fit one modifier at a time: age per 10 years (reference 60), Charlson per point (reference 2), complete six-hour modified SOFA per 2 points (reference 6), and 24-hour SOFA as a sensitivity. Sex compares male with female. Each recorded race category is compared separately with White; unknown/missing/unmapped values are excluded. Race is a recorded social classification, not a biological mechanism. The interaction OR is a ratio of pollution ORs, not a direct demographic or severity effect. Clinical and demographic interaction families have separate BH corrections retaining failed attempts. Six-hour SOFA is early post-entry organ dysfunction. Its interactions describe events of different observed severity and do not establish causal modification by baseline severity.</p>',
 '<p class="note">New modifier interactions currently use patient-cluster uncertainty. The baseline simulation calibration does not validate these interaction models. Complete-score and demographic exclusions may select the patients contributing to each contrast.</p>',
 '<h2>Analysis families</h2>',table(summary),
 '<h2>Weather adjustment sensitivity</h2>',table(weather_sensitivity[,.(organism,pollutant,exposure_window,n_events,OR=odds_ratio,CI_low=ci_low,CI_high=ci_high,q=q_value)]),
 '<p>Compare exposure_window_qc.csv to distinguish weather-related sample loss from covariate adjustment. Weather sensitivities have separate FDR families and are not primary findings.</p>',
 '<p>BH correction retains the full attempted hypothesis count in each family. Seven-day overall associations are primary; seven-day poverty interactions, other SES interactions, and diagnosis-defined outcomes have separate families. The 3-, 14-, and 28-day windows are corrected jointly within their corresponding sensitivity families.</p>',
 '<h2>Exposure duration sensitivity</h2>',table(fread(file.path(root,'exposure_window_qc.csv'))),
 '<p>Cases and referent rows are identical across windows within each pollutant and adjustment specification. Every retained row has complete daily exposure across all windows; only the weather sensitivity additionally requires complete lag-1–7 weather; each matched set has one case, at least one referent, and exposure variation in every window. The 28-day window overlaps heavily between nearby referent dates and may reduce precision. Different windows are correlated; differences in significance are not formal tests of differences in effects.</p>',
 '<img src="duration_sensitivity_forest.png" alt="Exposure window estimates for the six most frequent organisms" style="width:100%;height:auto">',table(duration),audit_html,
 '<h2>Culture selection and clinical sensitivities</h2>',table(clinical),
 '<p>Companion outcomes use all eligible early ICU admissions, any early respiratory culture, and any named organism. They help assess whether the organism signal parallels admission or testing patterns; they do not remove selection bias. Clinical sensitivities use 24/72-hour culture windows, pulmonary/upper-airway categories, and no documented inpatient antibacterial administration before culture. Absence of inpatient documentation does not establish absence of outpatient or referring-hospital treatment.</p>',
 '<h2>Shared calendar exposure inference checks</h2>',shared_html,
 '<p>Conditional Poisson estimates should agree with conditional logistic point estimates. Quasi-Poisson and calendar HAC (Bartlett weights, 28-day lag, scores aggregated across ZIPs) provide separate uncertainty checks; they are corrected separately by inference method. Sparse outcomes can invalidate all asymptotic methods, including dispersion estimates.</p>',
 '<h2>Simulation calibration</h2>',calibration_html,
 '<p>Rejection rates under OR=1 assess false positives; under OR=1.5 they assess conditional power. Report Monte Carlo intervals, failed fits and coverage alongside rates. Shared-calendar shocks are a sensitivity scenario, not a calibrated model of actual infection outbreaks. These simulations do not validate selection mechanisms or guarantee population-level power.</p>',
 '<h2>ACS linkage</h2>',table(fread(file.path(root,'ses_linkage_qc.csv'))),
 '<p>Fixed 2013–2017 ACS five-year estimates. Primary modifier: poverty, per 10 percentage points, centered at 20%. Secondary: no high school diploma (age 25+), unemployment (civilian labor force), household income, and crowding (&gt;1 occupants/room). Exact ZIP-to-ZCTA linkage is a geographic proxy; missing values are not imputed. Raw margins of error are retained but sampling uncertainty is not propagated into these models.</p>',
 '<h2>Admission diagnoses</h2>',table(fread(file.path(root,'diagnosis_capture.csv'))),
 '<p>Primary diagnosis flagged present on admission. Pneumonia/aspiration: ICD10 J12–J18/J69 or ICD9 480–486/507. Obstructive airway: ICD10 J41–J46 or ICD9 491/492/493/496. Remaining respiratory codes and nonrespiratory diagnoses are separate groups. Missing, ambiguous, and unmapped diagnoses remain explicit. Discharge coding is not proof of what was known at admission; groups define secondary outcomes.</p>',
 '<h2>Six most frequent estimable organisms</h2>',table(effects),
 '<p>For overall rows OR is per 5 ug/m3 PM2.5 or 10 ppb ozone. For poverty interaction rows OR is the ratio of pollution ORs per 10-point increase in neighborhood poverty; it is not the direct effect of poverty.</p>',
 '<h2>Federated outputs</h2><p>site_estimates.csv contains aggregate coefficients, standard errors, covariances, counts, and model status; no patient IDs, residential ZIPs, or dates. Pooling checks protocol compatibility and duplicate sites. Fixed-effect estimates are primary; heterogeneity and random-effect estimates are descriptive. Shared-exposure checks and simulations are development diagnostics; sparse results and departures from nominal error rates require resolution before definitive inference.</p>',
 '</body></html>')
writeLines(html,file.path(root,'report.html'))
message('Report: ',file.path(root,'report.html'))
