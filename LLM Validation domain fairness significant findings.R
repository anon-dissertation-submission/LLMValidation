# =============================================================================
# Equity and Fairness Assessment (Domain-Level)
# Assesses Fairness across Demographic Domains at the observation level.
# Outputs ONLY significant findings (p_fdr < 0.05) with descriptive domain names.
# =============================================================================

library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(purrr)
library(tools)

# Setup & Imports
source("bootstrap_helpers.R")

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"

N_BOOT     <- 10000  # 10000 iterations for stable p-values
SEED       <- 13

# Domain Mapping
# Updated to match the "Q1", "Q2" string format
domain_map <- tibble::tribble(
  ~q_num, ~domain_lab,
  "Q1",  "Unusual thoughts and experiences",
  "Q2",  "Suspiciousness",
  "Q3",  "Unusual somatic ideas",
  "Q4",  "Ideas of guilt",
  "Q5",  "Jealous ideas",
  "Q6",  "Unusual religious ideas",
  "Q7",  "Erotomanic ideas",
  "Q8",  "Grandiosity",
  "Q9",  "Auditory perceptual abnormalities",
  "Q10", "Visual perceptual abnormalities",
  "Q11", "Olfactory perceptual abnormalities",
  "Q12", "Gustatory perceptual abnormalities",
  "Q13", "Tactile perceptual abnormalities",
  "Q14", "Somatic perceptual abnormalities",
  "Q15", "Disorganised Communication Expression"
)

# Load and Prepare Domain-Level Data
cat("Loading and preparing domain-level data...\n")

# CAARMS thresholds for domain-level classification
CAARMS_SEVERITY_LOW_APS2a   <- 3
CAARMS_SEVERITY_HIGH_APS2a  <- 5
CAARMS_FREQUENCY_LOW_APS2a  <- 3
CAARMS_FREQUENCY_HIGH_APS2a <- 6
CAARMS_SEVERITY_APS2b       <- 6
CAARMS_FREQUENCY_APS2b      <- 3
CAARMS_SEVERITY_BLIPS       <- 6
CAARMS_FREQUENCY_LOW_BLIPS  <- 4
CAARMS_FREQUENCY_HIGH_BLIPS <- 6

is_chr_domain <- function(sev, freq) {
  (sev %in% CAARMS_SEVERITY_LOW_APS2a:CAARMS_SEVERITY_HIGH_APS2a &
     freq %in% CAARMS_FREQUENCY_LOW_APS2a:CAARMS_FREQUENCY_HIGH_APS2a) |
    (sev == CAARMS_SEVERITY_APS2b  & freq == CAARMS_FREQUENCY_APS2b) |
    (sev == CAARMS_SEVERITY_BLIPS  & freq %in% CAARMS_FREQUENCY_LOW_BLIPS:CAARMS_FREQUENCY_HIGH_BLIPS)
}

raw_data <- read_excel(DATA_FILE, sheet = DATA_SHEET) %>%
  rename(q_num = question) %>% # Renames the 'question' column to 'q_num'
  mutate(
    src_subject_id = as.character(src_subject_id),
    q_num = as.character(trimws(q_num)) # Ensures it is text and strips any hidden spaces
  )

# Calculate labels and join the domain_map
obs_df <- raw_data %>%
  filter(!is.na(severity), !is.na(frequency), !is.na(severity_pred), !is.na(frequency_pred), !is.na(q_num)) %>%
  left_join(domain_map, by = "q_num") %>%
  filter(!is.na(domain_lab)) %>% # Exclude any unmapped questions
  mutate(
    chr      = if_else(is_chr_domain(severity, frequency), "CHR", "HC"),
    chr_pred = if_else(is_chr_domain(severity_pred, frequency_pred), "CHR", "HC")
  )

# ------------------------------------------------------------------------------
# DEMOGRAPHICS: Safely check for columns without crashing
# ------------------------------------------------------------------------------
target_cols <- c("sex", "age_band", "language", "race", "site", "Site", "Race", "Sex", "Age")

demographics <- raw_data %>%
  select(src_subject_id, any_of(target_cols)) %>%
  distinct(src_subject_id, .keep_all = TRUE)

# Standardise names if they came in capitalised
if("Site" %in% names(demographics) && !"site" %in% names(demographics)) demographics <- rename(demographics, site = Site)
if("Race" %in% names(demographics) && !"race" %in% names(demographics)) demographics <- rename(demographics, race = Race)
if("Sex" %in% names(demographics) && !"sex" %in% names(demographics)) demographics <- rename(demographics, sex = Sex)

# Safely convert available columns to character to prevent factor errors
demographics <- demographics %>%
  mutate(across(any_of(c("sex", "age_band", "race", "site")), as.character))

# Process language grouping
if ("language" %in% names(demographics)) {
  demographics <- demographics %>%
    mutate(
      lang_group = case_when(
        is.na(language) ~ NA_character_,
        tolower(trimws(language)) == "english" ~ "English",
        TRUE ~ "Non-English"
      )
    )
}

# Process site grouping
if ("site" %in% names(demographics)) {
  demographics <- demographics %>%
    mutate(
      site_group = case_when(
        is.na(site) ~ NA_character_,
        toupper(trimws(site)) == "ME" ~ "ME", 
        TRUE ~ "Non-ME"
      )
    )
}

# Join demographics to observation data
domain_df <- obs_df %>%
  select(-any_of(c("sex", "age_band", "language", "lang_group", "language_binary", "race", "site", "site_group"))) %>% 
  left_join(demographics, by = "src_subject_id")


# Core Fairness Domain Bootstrap Function
domain_fairness_bootstrap <- function(df, d_lab, group_var, ref_group, comp_group, n_iter = 1000, seed = 13) {
  set.seed(seed)
  
  # Filter by demographic groups and the specific domain string label
  df_sub <- df %>% 
    filter(.data[[group_var]] %in% c(ref_group, comp_group) & domain_lab == d_lab) %>%
    filter(!is.na(chr) & !is.na(chr_pred))
  
  n_ref  <- sum(df_sub[[group_var]] == ref_group, na.rm = TRUE)
  n_comp <- sum(df_sub[[group_var]] == comp_group, na.rm = TRUE)
  
  # Skip if there's insufficient data to compare
  if (n_ref < 5 || n_comp < 5) return(NULL)
  
  unique_subj <- unique(df_sub$src_subject_id)
  n_subj      <- length(unique_subj)
  
  boot_res <- tibble(
    ref_ba=numeric(n_iter), comp_ba=numeric(n_iter), diff_ba=numeric(n_iter),
    ref_tpr=numeric(n_iter), comp_tpr=numeric(n_iter), diff_tpr=numeric(n_iter),
    ref_fpr=numeric(n_iter), comp_fpr=numeric(n_iter), diff_fpr=numeric(n_iter),
    ref_icc=numeric(n_iter), comp_icc=numeric(n_iter), diff_icc=numeric(n_iter),
    ref_r=numeric(n_iter),   comp_r=numeric(n_iter),   diff_r=numeric(n_iter),
    ref_mae=numeric(n_iter), comp_mae=numeric(n_iter), diff_mae=numeric(n_iter)
  )
  
  for (i in seq_len(n_iter)) {
    sampled_ids <- sample(unique_subj, n_subj, replace = TRUE)
    
    # Maintain subject-level clustering while resampling
    boot_sample <- df_sub %>%
      inner_join(
        tibble(src_subject_id = sampled_ids, .rep = seq_along(sampled_ids)),
        by = "src_subject_id", relationship = "many-to-many"
      ) %>% select(-.rep)
    
    # Binary Classification Metrics
    metrics_bin <- boot_sample %>%
      group_by(.data[[group_var]]) %>%
      group_modify(~ calc_metrics(.x)) %>%
      ungroup()
    
    # Continuous Reliability Metrics (using Severity as the domain proxy)
    metrics_cont <- boot_sample %>%
      group_by(.data[[group_var]]) %>%
      group_modify(~ {
        m <- calc_continuous_metrics(.x, "severity", "severity_pred")
        tibble(icc = m["icc"], r = m["r"], mae = m["mae"])
      }) %>%
      ungroup()
    
    ref_b <- metrics_bin %>% filter(.data[[group_var]] == ref_group)
    comp_b <- metrics_bin %>% filter(.data[[group_var]] == comp_group)
    ref_c <- metrics_cont %>% filter(.data[[group_var]] == ref_group)
    comp_c <- metrics_cont %>% filter(.data[[group_var]] == comp_group)
    
    if (nrow(ref_b) == 1 && nrow(comp_b) == 1 && nrow(ref_c) == 1 && nrow(comp_c) == 1) {
      # BA, TPR, FPR
      boot_res$ref_ba[i] <- ref_b$bal_acc; boot_res$comp_ba[i] <- comp_b$bal_acc; boot_res$diff_ba[i] <- comp_b$bal_acc - ref_b$bal_acc
      boot_res$ref_tpr[i] <- ref_b$sens;   boot_res$comp_tpr[i] <- comp_b$sens;   boot_res$diff_tpr[i] <- comp_b$sens - ref_b$sens
      boot_res$ref_fpr[i] <- 1-ref_b$spec; boot_res$comp_fpr[i] <- 1-comp_b$spec; boot_res$diff_fpr[i] <- (1-comp_b$spec) - (1-ref_b$spec)
      
      # ICC, Pearson, MAE
      boot_res$ref_icc[i] <- ref_c$icc; boot_res$comp_icc[i] <- comp_c$icc; boot_res$diff_icc[i] <- comp_c$icc - ref_c$icc
      boot_res$ref_r[i]   <- ref_c$r;   boot_res$comp_r[i]   <- comp_c$r;   boot_res$diff_r[i]   <- comp_c$r - ref_c$r
      boot_res$ref_mae[i] <- ref_c$mae; boot_res$comp_mae[i] <- comp_c$mae; boot_res$diff_mae[i] <- comp_c$mae - ref_c$mae
    } else {
      boot_res[i, ] <- NA
    }
  }
  
  # Format each metric to pivot easily
  format_metric <- function(metric_name, ref_v, comp_v, diff_v) {
    v <- na.omit(diff_v)
    if(length(v) < 2) return(NULL)
    
    p_val <- min(2 * min(mean(v <= 0), mean(v >= 0)), 1)
    
    tibble(
      Demographic = toTitleCase(group_var),
      Comparison = paste0(comp_group, " vs ", ref_group),
      `Symptom Domain` = d_lab, # Inject the nice string label here
      Metric = metric_name,
      `n (ref)` = n_ref,
      `n (comp)` = n_comp,
      
      ref_mean = mean(ref_v, na.rm=T), ref_lo = quantile(ref_v, 0.025, na.rm=T, names=F), ref_hi = quantile(ref_v, 0.975, na.rm=T, names=F),
      comp_mean = mean(comp_v, na.rm=T), comp_lo = quantile(comp_v, 0.025, na.rm=T, names=F), comp_hi = quantile(comp_v, 0.975, na.rm=T, names=F),
      diff_mean = mean(v), diff_lo = quantile(v, 0.025, names=F), diff_hi = quantile(v, 0.975, names=F),
      p_val = p_val
    )
  }
  
  bind_rows(
    format_metric("BA", boot_res$ref_ba, boot_res$comp_ba, boot_res$diff_ba),
    format_metric("TPR", boot_res$ref_tpr, boot_res$comp_tpr, boot_res$diff_tpr),
    format_metric("FPR", boot_res$ref_fpr, boot_res$comp_fpr, boot_res$diff_fpr),
    format_metric("ICC", boot_res$ref_icc, boot_res$comp_icc, boot_res$diff_icc),
    format_metric("Pearson", boot_res$ref_r, boot_res$comp_r, boot_res$diff_r),
    format_metric("MAE", boot_res$ref_mae, boot_res$comp_mae, boot_res$diff_mae)
  )
}

# Setup Dynamic Comparisons
cat("\nConfiguring dynamic demographic comparisons across domains...\n")

get_comparisons <- function(df, var_name, forced_ref = NULL) {
  counts <- df %>% filter(!is.na(.data[[var_name]])) %>% count(.data[[var_name]], sort = TRUE)
  if (nrow(counts) < 2) return(list()) 
  
  ref <- if (!is.null(forced_ref) && forced_ref %in% counts[[var_name]]) forced_ref else counts[[var_name]][1] 
  comps <- setdiff(counts[[var_name]], ref)
  
  lapply(comps, function(c) list(var = var_name, ref = ref, comp = c))
}

comparisons_to_run <- list()
if ("sex" %in% names(domain_df))        comparisons_to_run <- c(comparisons_to_run, get_comparisons(domain_df, "sex", "Male"))
if ("age_band" %in% names(domain_df))   comparisons_to_run <- c(comparisons_to_run, get_comparisons(domain_df, "age_band", "18-25"))
if ("lang_group" %in% names(domain_df)) comparisons_to_run <- c(comparisons_to_run, get_comparisons(domain_df, "lang_group", "English"))
if ("race" %in% names(domain_df))       comparisons_to_run <- c(comparisons_to_run, get_comparisons(domain_df, "race"))
if ("site_group" %in% names(domain_df)) comparisons_to_run <- c(comparisons_to_run, get_comparisons(domain_df, "site_group", "ME"))

if (length(comparisons_to_run) == 0) stop("No valid demographic columns found!")

# Get all unique domain labels from the joined map
unique_domains <- unique(na.omit(domain_df$domain_lab))

# Run Bootstraps (Nested loop for Domain x Demographic)
cat(sprintf("Running domain-level comparisons (Total domains: %d). This will take a while...\n", length(unique_domains)))

results_raw_list <- list()
counter <- 1

for (domain in unique_domains) {
  for (cmp in comparisons_to_run) {
    cat(sprintf("\rEvaluating Domain: %-35s | %s: %s vs %s", domain, cmp$var, cmp$comp, cmp$ref))
    
    res <- domain_fairness_bootstrap(
      domain_df, d_lab = domain, group_var = cmp$var, 
      ref_group = cmp$ref, comp_group = cmp$comp, n_iter = N_BOOT, seed = SEED
    )
    
    if(!is.null(res)) {
      results_raw_list[[counter]] <- res
      counter <- counter + 1
    }
  }
}

cat("\n\nBootstrap complete! Formatting results...\n")
all_results <- bind_rows(results_raw_list)

# FDR Adjustment & Format Final Tables

# Group by Metric to apply FDR correction independently across the metric space
adjusted_results <- all_results %>%
  group_by(Metric) %>%
  mutate(p_fdr = p.adjust(p_val, method = "fdr")) %>%
  ungroup()

fmt_ci <- function(m, lo, hi) {
  if (is.na(m)) return("NA")
  sprintf("%.3f [%.3f, %.3f]", m, lo, hi)
}
fmt_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<.001")
  sprintf("%.3f", p)
}

formatted_table <- adjusted_results %>%
  rowwise() %>%
  mutate(
    `Ref Value [95% CI]`  = fmt_ci(ref_mean, ref_lo, ref_hi),
    `Comp Value [95% CI]` = fmt_ci(comp_mean, comp_lo, comp_hi),
    `Δ [95% CI]`          = fmt_ci(diff_mean, diff_lo, diff_hi),
    p_fdr                 = fmt_pval(p_fdr)
  ) %>%
  ungroup() %>%
  select(
    Demographic, Comparison, `Symptom Domain`, Metric,
    `n (ref)`, `n (comp)`,
    `Ref Value [95% CI]`, `Comp Value [95% CI]`, `Δ [95% CI]`, p_fdr
  )

# Sort by Demographic, Comparison, Symptom Domain, then Metric
final_table <- formatted_table %>%
  arrange(Demographic, Comparison, `Symptom Domain`, Metric)

# Save the full p-values background data just in case
write_csv(adjusted_results, "all_pvalues_domain.csv")

# Filter ONLY significant findings (p_fdr < 0.05) and save
significant_findings <- final_table %>%
  filter(!is.na(p_fdr) & p_fdr != "NA") %>%
  filter(as.numeric(gsub("<", "", p_fdr)) < 0.05)

write_csv(significant_findings, "domain_fairness_significant_findings.csv")

cat("\nDone! Saved all combinations to 'all_pvalues_domain.csv'\n")
cat(sprintf("Saved %d significant findings to 'domain_fairness_significant_findings.csv'\n\n", nrow(significant_findings)))

print(knitr::kable(head(significant_findings, 10)))