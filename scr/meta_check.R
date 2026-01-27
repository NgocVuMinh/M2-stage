
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# Initialize result dataframes
matched_df <- data.frame()
mismatch_tissue <- data.frame()
mismatch_stage <- data.frame()
mismatch_sampling_site <- data.frame()
mismatch_status <- data.frame()

err_ids <- sel$Run

# Loop through each ERR ID
for (err_id in err_ids) {
  
  # Get metadata from sel for this ERR ID
  ebi_row <- sel[sel$Run == err_id, ]
  
  # Extract patient ID and metadata from EBI
  patient_id <- ebi_row$Sample.Characteristic.individual.[1]
  ebi_tissue <- ebi_row$Sample.Characteristic.organism.part.[1]
  ebi_stage <- ebi_row$pcw_stage[1]
  ebi_sampling_site <- ebi_row$Sample.Characteristic.sampling.site.[1]
  
  # Find all rows in sel.supp for this patient
  supp_rows <- sel.supp[sel.supp$id_hdbr == patient_id, ]
  
  # Skip if patient not found in sel.supp
  if (nrow(supp_rows) == 0) {
    print(paste0("Patient not found in supp: ", patient_id))
    next
  }
  
  # Initialize mismatch flags
  tissue_match <- FALSE
  stage_match <- FALSE
  sampling_site_match <- FALSE
  
  # Check for perfect matches (all 3 factors)
  perfect_matches <- supp_rows %>%
    filter(tissue_clean == ebi_tissue,
           pcw_stage == ebi_stage,
           left_right == ebi_sampling_site)
  
  # If perfect match exists (one or more), save to matched_df
  if (nrow(perfect_matches) > 0) {
    matched_df <- rbind(matched_df, ebi_row)
    tissue_match <- TRUE
    stage_match <- TRUE
    sampling_site_match <- TRUE
  } else {
    # Check partial matches to identify mismatches
    
    # Check tissue matches
    tissue_candidates <- supp_rows %>%
      filter(tissue_clean == ebi_tissue)
    
    if (nrow(tissue_candidates) > 0) {
      tissue_match <- TRUE
      
      # Among tissue matches, check stage
      stage_candidates <- tissue_candidates %>%
        filter(pcw_stage == ebi_stage)
      
      if (nrow(stage_candidates) > 0) {
        stage_match <- TRUE
        
        # Among stage matches, check sampling site
        site_candidates <- stage_candidates %>%
          filter(left_right == ebi_sampling_site)
        
        if (nrow(site_candidates) > 0) {
          sampling_site_match <- TRUE
        } else {
          # Sampling site mismatch
          mismatch_sampling_site <- rbind(mismatch_sampling_site,
                                          data.frame(
                                            ERR_ID = err_id,
                                            patient_ID = patient_id,
                                            sampling_site_EBI = ebi_sampling_site,
                                            sampling_site_SUPP = paste(unique(stage_candidates$left_right), collapse = "; ")
                                          ))
        }
      } else {
        # Stage mismatch
        mismatch_stage <- rbind(mismatch_stage,
                                data.frame(
                                  ERR_ID = err_id,
                                  patient_ID = patient_id,
                                  stage_EBI = ebi_stage,
                                  stage_SUPP = paste(unique(tissue_candidates$pcw_stage), collapse = "; ")
                                ))
      }
    } else {
      # Tissue mismatch
      mismatch_tissue <- rbind(mismatch_tissue,
                               data.frame(
                                 ERR_ID = err_id,
                                 patient_ID = patient_id,
                                 tissue_EBI = ebi_tissue,
                                 tissue_SUPP = paste(unique(supp_rows$tissue_clean), collapse = "; ")
                               ))
    }
  }
  
  # Record overall mismatch status
  mismatch_status <- rbind(mismatch_status,
                           data.frame(
                             ERR_ID = err_id,
                             patient_ID = patient_id,
                             tissue = ifelse(tissue_match, 0, 1),
                             stage = ifelse(stage_match, 0, 1),
                             sampling_site = ifelse(sampling_site_match, 0, 1)
                           ))
}

# View results
print(paste("Total samples analyzed:", nrow(mismatch_status)))
print(paste("Perfect matches:", nrow(matched_df)))
print(paste("Tissue mismatches:", nrow(mismatch_tissue)))
print(paste("Stage mismatches:", nrow(mismatch_stage)))
print(paste("Sampling site mismatches:", nrow(mismatch_sampling_site)))

# View mismatch summaries
head(mismatch_tissue)
head(mismatch_stage)
head(mismatch_sampling_site)
head(mismatch_status)