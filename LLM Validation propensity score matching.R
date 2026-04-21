# Suppress warnings to keep output clean
options(warn = -1)
# =============================================================================
# Propensity Score Matching (PSM) Fairness Assessment
# Method: Nearest-Neighbour Matching (Ho et al., 2011)
# =============================================================================

library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(MatchIt)
library(purrr)
library(future)
library(furrr)
library(ggplot2)


# Setup & Imports
source("bootstrap_helpers.R")

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"
N_BOOT     <- 10 # Updated to 10 for testing
SEED       <- 13

# Classification Logic
is_chr_domain <- function(sev, freq) {
  (sev %in% 3:5 & freq %in% 3:6) | # APS2a
    (sev == 6 & freq == 3) |         # APS2b
    (sev == 6 & freq %in% 4:6)       # BLIPS
}

# Data Preparation (Visit-Level)
cat("Preparing visit-level data with balancing covariates...\n")

raw_data <- read_excel(DATA_FILE, sheet = DATA_SHEET) %>%
  rename_with(~tolower(.), any_of(c("Sex", "Race", "Site", "Age", "Language", "age_band"))) %>%
  mutate(src_subject_id = as.character(src_subject_id))

# Create observation-level labels
obs_df <- raw_data %>%
  filter(!is.na(severity), !is.na(frequency), !is.na(severity_pred), !is.na(frequency_pred)) %>%
  mutate(
    chr      = if_else(is_chr_domain(severity, frequency), "CHR", "HC"),
    chr_pred = if_else(is_chr_domain(severity_pred, frequency_pred), "CHR", "HC")
  )

# Aggregate to visit level with clinical & demographic covariates
visit_df <- obs_df %>%
  group_by(src_subject_id, visit) %>%
  summarize(
    chr = if_else(any(chr == "CHR"), "CHR", "HC"),
    chr_pred = if_else(any(chr_pred == "CHR"), "CHR", "HC"),
    
    mean_severity = mean(severity, na.rm = TRUE),
    mean_frequency = mean(frequency, na.rm = TRUE),
    
    sex = first(sex),
    race = first(race),
    site = first(site),
    age_band = first(age_band),
    language = first(language),
    .groups = "drop"
  ) %>%
  mutate(
    chr_numeric = if_else(chr == "CHR", 1, 0),
    # Updated logic: Separate 'Unavailable' from 'Non-English'
    language_bin = case_when(
      tolower(trimws(language)) == "english" ~ "English",
      tolower(trimws(language)) == "unavailable" | is.na(language) ~ "Unavailable",
      TRUE ~ "Non-English"
    )
  )

# Core PSM Function (Upgraded)
run_psm_analysis <- function(df, group_var, ref_val, comp_val, match_vars) {
  
  # Standardize case to prevent text mismatches
  ref_val_clean <- tolower(trimws(ref_val))
  comp_val_clean <- tolower(trimws(comp_val))
  
  # Filter, Drop NAs, and Assign Binary Treatment
  df_sub <- df %>% 
    mutate(temp_group = tolower(trimws(as.character(.data[[group_var]])))) %>%
    filter(temp_group %in% c(ref_val_clean, comp_val_clean)) %>%
    drop_na(any_of(match_vars)) %>%
    mutate(treated = as.numeric(temp_group == comp_val_clean))
  
  # Failsafes to prevent matchit() crashes
  if (nrow(df_sub) < 10) {
    cat(sprintf("\n[Skipped] %s vs %s: Not enough total participants with complete data.\n", comp_val, ref_val))
    return(NULL)
  }
  
  if (length(unique(df_sub$treated)) != 2) {
    cat(sprintf("\n[Skipped] %s vs %s: One of the groups has 0 participants after dropping missing covariates.\n", comp_val, ref_val))
    return(NULL)
  }
  
  # A) UNMATCHED METRICS
  unmatched_res <- binary_bootstrap(df_sub, "temp_group", comp_val_clean, ref_val_clean, 
                                    unique(df_sub$src_subject_id), 
                                    length(unique(df_sub$src_subject_id)), N_BOOT)
  
  # B) NEAREST-NEIGHBOUR MATCHING WITH CALIPER (Ho et al., 2011)
  psm_formula <- as.formula(paste("treated ~", paste(match_vars, collapse = " + ")))
  
  # Added a 0.2 SD caliper to enforce strict similarity between matched pairs
  match_obj <- matchit(psm_formula, data = df_sub, 
                       method = "nearest", 
                       distance = "glm", 
                       replace = FALSE, 
                       ratio = 1,
                       caliper = 0.2) 
  
  df_matched <- match.data(match_obj)
  
  # Ensure we still have enough matched pairs after applying the strict caliper
  if (nrow(df_matched) < 10) {
    cat(sprintf("\n[Skipped] %s vs %s: Not enough valid matches found within the 0.2 caliper limit.\n", comp_val, ref_val))
    return(NULL)
  }
  
  # C) AUTOMATED COVARIATE BALANCE CHECK
  bal_summary <- summary(match_obj, un = FALSE)$sum.matched
  
  # Locate the Standardized Mean Difference column
  smd_col <- grep("Std\\. Mean Diff", colnames(bal_summary), value = TRUE)
  max_smd <- NA
  
  if (length(smd_col) > 0) {
    # Isolate covariates (exclude the underlying propensity score distance metric)
    cov_smds <- bal_summary[rownames(bal_summary) != "distance", smd_col]
    max_smd <- max(abs(cov_smds), na.rm = TRUE)
  }
  
  # Flag any matches where the maximum SMD exceeds the standard 0.1 threshold
  balance_flag <- if_else(!is.na(max_smd) && max_smd > 0.1, "WARNING: Poor Balance (>0.1)", "Balanced")
  
  # D) MATCHED METRICS
  matched_res <- binary_bootstrap(df_matched, "temp_group", comp_val_clean, ref_val_clean, 
                                  unique(df_matched$src_subject_id), 
                                  length(unique(df_matched$src_subject_id)), N_BOOT)
  
  # E) FORMAT OUTPUT
  get_row <- function(metric_key, metric_label) {
    u <- unmatched_res$summaries[[metric_key]]
    m <- matched_res$summaries[[metric_key]]
    
    u_sig <- !(u$lo <= 0 & u$hi >= 0)
    m_sig <- !(m$lo <= 0 & m$hi >= 0)
    
    interp <- case_when(
      !u_sig & !m_sig ~ "Not significant in either",
      u_sig & !m_sig  ~ "Attenuated after matching",
      u_sig & m_sig   ~ "Persists after matching",
      TRUE            ~ "NS after matching"
    )
    
    tibble(
      Comparison = paste(comp_val, "vs", ref_val),
      `n pairs`  = nrow(df_matched) / 2,
      Metric     = paste("Δ", metric_label),
      `Unmatched Δ [95% CI]` = sprintf("%.3f [%.3f, %.3f]", u$mean, u$lo, u$hi),
      `Matched Δ [95% CI]`   = sprintf("%.3f [%.3f, %.3f]", m$mean, m$lo, m$hi),
      Interpretation = interp,
      `Max SMD` = if (!is.na(max_smd)) sprintf("%.3f", max_smd) else "N/A",
      `Balance Status` = balance_flag
    )
  }
  
  bind_rows(
    get_row("bal_acc", "BA"),
    get_row("sens", "TPR"),
    get_row("spec", "FPR")
  )
}

# Execution
cat("\nRunning Matching Analysis (Nearest-Neighbour) in PARALLEL...\n")

common_vars <- c("site", "age_band", "sex", "mean_severity", "mean_frequency")

tasks <- list(
  list(var = "race", ref = "White", comp = "Black or African American", match = common_vars),
  list(var = "race", ref = "White", comp = "Asian", match = common_vars),
  list(var = "age_band", ref = "18-25", comp = "<18", match = setdiff(common_vars, "age_band")),
  list(var = "language_bin", ref = "English", comp = "Non-English", match = common_vars),
  list(var = "language_bin", ref = "English", comp = "Unavailable", match = common_vars)
)

plan(multisession, workers = availableCores())

# Added .progress = TRUE to generate a live progress bar in the console
results <- future_map(tasks, function(t) {
  run_psm_analysis(visit_df, t$var, t$ref, t$comp, t$match)
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Filter out any NULL results from skipped comparisons
results <- compact(results)

if (length(results) > 0) {
  final_table <- bind_rows(results) %>%
    mutate(Metric = if_else(Metric == "Δ SPEC", "Δ FPR", Metric))
  
  write_csv(final_table, "psm_matched_fairness_metrics.csv")
  cat("\nAnalysis complete. Results saved to 'psm_matched_fairness_metrics.csv'.\n")
  print(knitr::kable(final_table))
  
}