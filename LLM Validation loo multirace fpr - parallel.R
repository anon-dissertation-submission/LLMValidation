library(dplyr)
library(purrr)
library(readr)
library(readxl) 
library(irr)
library(parallel) # Base R parallelization
library(pbapply)  # Bulletproof progress bars for parallel clusters

source("bootstrap_helpers.R")
source("visit_level_helpers.R")

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"
N_BOOTS <- 10000  # Keep this at 10 for your first test!

cat("Loading and preprocessing data...\n")
phenotype_data <- read_excel(DATA_FILE, sheet = DATA_SHEET)

phenotype_data <- phenotype_data %>%
  mutate(
    sex            = factor(sex),
    race           = factor(race),
    site           = factor(substr(src_subject_id, 1, 2)),
    src_subject_id = as.character(src_subject_id),
    age_band       = factor(age_band),
    language       = factor(language),
    language_binary = case_when(
      tolower(trimws(language)) == "english" ~ "English",
      is.na(language) ~ NA_character_,
      TRUE ~ "Non-English"
    )
  )

visit_data <- collapse_to_visit_level(phenotype_data)

run_bias_bootstrap <- function(data, group_var, ref_group = NULL, n_iterations = 1000, seed = 13) {
  set.seed(seed)
  if (!group_var %in% names(data)) stop(paste("Column", group_var, "not found in data."))
  
  unique_groups   <- unique(as.character(data[[group_var]]))
  unique_groups   <- unique_groups[!is.na(unique_groups)] 
  
  unique_subjects <- unique(data$src_subject_id)
  unique_subjects <- unique_subjects[!is.na(unique_subjects)]
  n_subjects      <- length(unique_subjects)
  
  is_binary_mode <- (length(unique_groups) == 2 && is.null(ref_group))
  
  if (is_binary_mode) {
    g1 <- unique_groups[1]
    g2 <- unique_groups[2]
    return(binary_bootstrap(data, group_var, g1, g2, unique_subjects, n_subjects, n_iterations))
  } else {
    comp_groups <- setdiff(unique_groups, ref_group)
    return(multigroup_bootstrap(data, group_var, ref_group, comp_groups, unique_subjects, n_subjects, n_iterations))
  }
}

multirace_hc_subjects <- visit_data %>%
  filter(race == "More than one race", chr == "HC") %>%
  pull(src_subject_id) %>%
  unique()

cat(sprintf("Running LOO for %d Multi-race HC subjects...\n", length(multirace_hc_subjects)))

# =============================================================================
# Windows Parallel Cluster Setup (The pbapply method)
# =============================================================================
available_workers <- max(1, parallel::detectCores() - 1)

cat("\nSpinning up Windows background workers... (If it hangs here, check your Firewall)\n")
cl <- makeCluster(available_workers)

# Export variables from your main session to the background workers
clusterExport(cl, varlist = c("visit_data", "run_bias_bootstrap", "N_BOOTS", "specdiff_to_fprdiff"))

# Force the background workers to load packages and source helper files
clusterEvalQ(cl, {
  library(dplyr)
  library(irr)
  library(tidyr)
  library(rlang)
  
  # Source these directly into the worker's brain
  source("bootstrap_helpers.R", local = TRUE)
  source("visit_level_helpers.R", local = TRUE)
})

cat("Workers ready! Starting parallel execution...\n")

# pbapply replaces purrr and automatically handles a clean progress bar
# It also safely wraps the random seed generation for parallel clusters
results_list <- pbapply::pblapply(multirace_hc_subjects, function(dropped_subj) {
  
  dat <- visit_data %>% filter(src_subject_id != dropped_subj)
  res <- run_bias_bootstrap(dat, group_var = "race", ref_group = "White",
                            n_iterations = N_BOOTS)
  
  fpr_df <- specdiff_to_fprdiff(res$summaries$spec) %>%
    filter(group == "More than one race") %>%
    mutate(dropped_subject = dropped_subj)
  
  return(fpr_df)
  
}, cl = cl) 

# SHUT DOWN THE CLUSTER (Crucial step so you don't leak RAM)
stopCluster(cl)

# Combine the list of results back into a dataframe
loo_multirace <- bind_rows(results_list)

# =============================================================================
# Baseline Calculation (Sequential)
# =============================================================================
cat("\nRunning baseline (full sample)...\n")
res_baseline <- run_bias_bootstrap(visit_data, group_var = "race", 
                                   ref_group = "White", n_iterations = N_BOOTS)

baseline <- specdiff_to_fprdiff(res_baseline$summaries$spec) %>%
  filter(group == "More than one race") %>%
  mutate(dropped_subject = "None (full sample)")

loo_full <- bind_rows(baseline, loo_multirace)
write_csv(loo_full, "loo_multirace_fpr.csv")

cat("\nComplete! Results saved to loo_multirace_fpr.csv\n")
print(loo_full)
