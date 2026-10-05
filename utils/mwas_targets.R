suppressPackageStartupMessages(library(data.table))
mwas_target_ids <- function(b,organism,analysis='overall') {
 window <- if(analysis=='clinical:culture24')24 else if(analysis=='clinical:culture72')72 else 48
 e <- b$episodes[hours_from_icu<=window];det <- b$detections[hours_from_icu<=window]
 if(organism=='__early_icu_admission')return(b$cohort$hospitalization_id)
 if(organism=='__any_respiratory_culture')return(unique(e$hospitalization_id))
 if(organism=='__any_named_organism')return(unique(det$hospitalization_id))
 if(analysis=='clinical:no_documented_antibiotics')det <- det[ab_status=='none_documented']
 if(analysis %in% c('clinical:pulmonary_sources','clinical:upper_airway_sources')) {
  sources <- if(analysis=='clinical:pulmonary_sources')c('respiratory_tract','respiratory_tract_lower') else
    c('nasopharynx_upperairway','oropharynx_tongue_oralcavity')
  det <- det[episode_id %in% e[fluid_category %in% sources,episode_id]]
 }
 unique(det[organism_category==organism,hospitalization_id])
}
