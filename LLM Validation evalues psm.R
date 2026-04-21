# Install required packages if needed
# install.packages(c("readr", "dplyr", "stringr", "tidyr"))

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# Helper function to compute the E-value from a Risk Ratio (RR)
calc_evalue <- function(rr) {
  if (is.na(rr)) return(NA)
  if (rr < 1) rr <- 1 / rr 
  if (rr == 1) return(1)
  
  return(rr + sqrt(rr * (rr - 1)))
}
calc_evalue_vec <- Vectorize(calc_evalue)

# Load the datasets
df_unmatched <- read_csv("visit_fairness_summary_10000Boots_331.csv", show_col_types = FALSE)
df_psm <- read_csv("psm_matched_fairness_metrics_fdr.csv", show_col_types = FALSE)

# Helper function to parse "Estimate [LCL, UCL]" strings
parse_ci_string <- function(df, col_name, prefix) {
  df %>%
    extract(
      col = !!sym(col_name),
      into = c(paste0(prefix, "_Est"), paste0(prefix, "_LCL"), paste0(prefix, "_UCL")),
      regex = "([\\-\\.\\d]+)\\s*\\[\\s*([\\-\\.\\d]+),\\s*([\\-\\.\\d]+)\\s*\\]",
      remove = FALSE,
      convert = TRUE
    )
}

# Parse unmatched (baseline) risks
df_unmatched_parsed <- df_unmatched %>%
  parse_ci_string("Ref TPR [95% CI]", "Ref_TPR") %>%
  parse_ci_string("Ref FPR [95% CI]", "Ref_FPR") %>%
  select(Comparison, Ref_TPR_Est, Ref_FPR_Est)

# Parse matched differences (PSM)
df_psm_parsed <- df_psm %>%
  parse_ci_string("Matched Δ [95% CI]", "Matched_Diff") %>%
  filter(Metric %in% c("Δ TPR", "Δ FPR"))

# Map "vs" format (Matched) to the baseline Category/Reference format (Unmatched)
df_psm_mapped <- df_psm_parsed %>%
  mutate(
    Unmatched_Comp = case_when(
      grepl("Black", Comparison) ~ "Race: Black or African American - White",
      grepl("Asian", Comparison) ~ "Race: Asian - White",
      grepl("<18", Comparison) ~ "Age_band: <18 - 18-25",
      grepl("Non-English", Comparison) ~ "Lang_group: Non-English - English",
      
      # FIX: Map "Unavailable vs English" to the row containing the "English" reference risks
      grepl("Unavailable", Comparison) ~ "Lang_group: Non-English - English", 
      
      TRUE ~ NA_character_
    )
  )

# Join datasets, compute exact Matched RRs, and calculate E-values
evalue_results <- df_psm_mapped %>%
  left_join(df_unmatched_parsed, by = c("Unmatched_Comp" = "Comparison")) %>%
  mutate(
    # Select the proper baseline risk based on metric
    Ref_Risk = if_else(Metric == "Δ TPR", Ref_TPR_Est, Ref_FPR_Est),
    
    # Calculate the Matched Risk Ratio
    Comp_Risk = pmax(0.0001, pmin(0.9999, Ref_Risk + Matched_Diff_Est)),
    Matched_RR = Comp_Risk / Ref_Risk,
    
    # Calculate E-values
    E_Value_Point = calc_evalue_vec(Matched_RR),
    
    # Assess if CI crosses the null
    Significant = !(Matched_Diff_LCL <= 0 & Matched_Diff_UCL >= 0),
    
    # If the CI crosses zero, the boundary E-value is strictly 1.0
    E_Value_CI = if_else(
      !Significant, 
      1.0, 
      pmin(
        calc_evalue_vec((Ref_Risk + Matched_Diff_LCL) / Ref_Risk),
        calc_evalue_vec((Ref_Risk + Matched_Diff_UCL) / Ref_Risk)
      )
    )
  ) %>%
  select(
    `PSM Comparison` = Comparison,
    Metric,
    `Matched Δ [95% CI]`,
    Ref_Risk,
    Matched_RR,
    E_Value_Point,
    E_Value_CI
  )

# View the final, corrected results
print(evalue_results, width = Inf)
write_csv(evalue_results, "evalues_psm.csv")
