# ==============================================================================
# UNIFIED MASTER SCRIPT: PROTEOME-WIDE MENDELIAN RANDOMIZATION 
# ==============================================================================
# Description: A scalable function that automates TwoSampleMR, pleiotropy checks, 
# FDR correction, and Excel reporting for any input cardiovascular phenotype.
# ==============================================================================

library(TwoSampleMR)
library(openxlsx)
library(dplyr)

# ==============================================================================
# 1. DEFINE THE CORE PIPELINE FUNCTION
# ==============================================================================
run_szokodi_pipeline <- function(harmonized_dat, outcome_name) {
  cat("\n===================================================\n")
  cat("INITIATING PIPELINE FOR:", outcome_name, "\n")
  cat("===================================================\n")
  
  # Ensure ID tags are properly mapped
  harmonized_dat$id.exposure <- harmonized_dat$exposure
  harmonized_dat$id.outcome <- outcome_name
  
  # Main MR & Odds Ratios
  cat("Running MR math and calculating Odds Ratios...\n")
  mr_res <- mr(harmonized_dat)
  mr_res <- generate_odds_ratios(mr_res)
  mr_res$P.adjust <- p.adjust(mr_res$pval, method = "BH")
  
  # Sensitivity Analyses
  cat("Running Sensitivity Diagnostics...\n")
  het_res <- mr_heterogeneity(harmonized_dat)
  pleio_res <- mr_pleiotropy_test(harmonized_dat)
  
  # Format Sheets
  cat("Formatting Szokodi Excel template...\n")
  sheet1_mr <- mr_res %>% 
    select(Exposure = exposure, id.outcome, Method = method, nsnp, b, 
           lo_ci, up_ci, se, `P-val` = pval, OR = or, OR_lci95 = or_lci95, 
           OR_uci95 = or_uci95, P.adjust) %>%
    arrange(P.adjust)
    
  sheet2_het <- if (nrow(het_res) > 0) {
    het_res %>% select(Exposure = exposure, Method = method, Q, Q_df, Q_pval)
  } else { data.frame(Message = "No heterogeneity data") }
    
  sheet3_pleio <- if (nrow(pleio_res) > 0) {
    pleio_res %>% select(id.exposure, id.outcome, Outcome = outcome, Exposure = exposure, 
                         Egger_intercept = egger_intercept, se, pval)
  } else { data.frame(Message = "No pleiotropy data") }
  
  # Export multi-sheet workbook
  if (!dir.exists("results")) dir.create("results")
  wb <- createWorkbook()
  addWorksheet(wb, "Confirmatory MR")
  writeData(wb, "Confirmatory MR", sheet1_mr)
  addWorksheet(wb, "Heterogeneity")
  writeData(wb, "Heterogeneity", sheet2_het)
  addWorksheet(wb, "Pleiotropy")
  writeData(wb, "Pleiotropy", sheet3_pleio)
  
  save_path <- paste0("results/Results_for_", outcome_name, ".xlsx")
  saveWorkbook(wb, save_path, overwrite = TRUE)
  
  cat("SUCCESS! Report saved to:", save_path, "\n")
}

# ==============================================================================
# 2. EXECUTE PIPELINE ACROSS MULTIPLE PHENOTYPES
# ==============================================================================
# Note: This assumes dat_hcm_decode and dat_dcm_decode have already been 
# extracted and harmonized via the standard TwoSampleMR process in your environment.

# Run Hypertrophic Cardiomyopathy
run_szokodi_pipeline(harmonized_dat = dat_hcm_decode, outcome_name = "HCM_deCODE")

# Run Dilated Cardiomyopathy
run_szokodi_pipeline(harmonized_dat = dat_dcm_decode, outcome_name = "DCM_deCODE")