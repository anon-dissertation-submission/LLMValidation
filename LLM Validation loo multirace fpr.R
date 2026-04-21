# Leave-One-Out (LOO) cross-validation analyzing false positive rates (FPR) across racial groups.
# Assesses how model generalizes when trained on all races except one, testing robustness across diverse populations.
library(dplyr)
library(purrr)
library(readr)
library(irr)
source("bootstrap_helpers.R")
source("visit_level_helpers.R")

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"
N_BOOTS <- 10000

phenotype_data <- read_excel(DATA_FILE, sheet = DATA_SHEET)

# Apply standard preprocessing
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

# Router: binary vs multigroup bootstrap
run_bias_bootstrap <- function(data, group_var, ref_group = NULL, n_iterations = 1000, seed = 13) {
  set.seed(seed)
  
  if (!group_var %in% names(data)) stop(paste("Column", group_var, "not found in data."))
  
  # STRIP NA VALUES HERE SO THEY NEVER BECOME A "GROUP"
  unique_groups   <- unique(as.character(data[[group_var]]))
  unique_groups   <- unique_groups[!is.na(unique_groups)] 
  
  unique_subjects <- unique(data$src_subject_id)
  unique_subjects <- unique_subjects[!is.na(unique_subjects)]
  n_subjects      <- length(unique_subjects)
  
  is_binary_mode <- (length(unique_groups) == 2 && is.null(ref_group))
  
  if (is_binary_mode) {
    g1 <- unique_groups[1]
    g2 <- unique_groups[2]
    cat(sprintf("Starting BINARY bootstrap for %s (%s - %s) with %d iterations...\n",
                toupper(group_var), g1, g2, n_iterations))
    return(binary_bootstrap(data, group_var, g1, g2, unique_subjects, n_subjects, n_iterations))
  } else {
    if (is.null(ref_group)) stop("You must provide a 'ref_group' for variables with more than 2 groups.")
    if (!ref_group %in% unique_groups) stop(paste("Reference group", ref_group, "not found."))
    
    comp_groups <- setdiff(unique_groups, ref_group)
    cat(sprintf("Starting MULTI-GROUP bootstrap for %s (Ref: %s) with %d iterations...\n",
                toupper(group_var), ref_group, n_iterations))
    return(multigroup_bootstrap(data, group_var, ref_group, comp_groups, unique_subjects, n_subjects, n_iterations))
  }
}

multirace_hc_subjects <- visit_data %>%
  filter(race == "More than one race", chr == "HC") %>%
  pull(src_subject_id) %>%
  unique()

cat(sprintf("Running LOO for %d Multi-race HC subjects\n", length(multirace_hc_subjects)))

loo_multirace <- map_dfr(multirace_hc_subjects, function(dropped_subj) {
  dat <- visit_data %>% filter(src_subject_id != dropped_subj)
  res <- run_bias_bootstrap(dat, group_var = "race", ref_group = "White",
                            n_iterations = N_BOOTS)
  fpr_df <- specdiff_to_fprdiff(res$summaries$spec) %>%
    filter(group == "More than one race") %>%
    mutate(dropped_subject = dropped_subj)
  return(fpr_df)
})

# Also run with no subjects dropped for comparison baseline
res_baseline <- run_bias_bootstrap(visit_data, group_var = "race", 
                                   ref_group = "White", n_iterations = N_BOOTS)
baseline <- specdiff_to_fprdiff(res_baseline$summaries$spec) %>%
  filter(group == "More than one race") %>%
  mutate(dropped_subject = "None (full sample)")

loo_full <- bind_rows(baseline, loo_multirace)
write_csv(loo_full, "loo_multirace_fpr.csv")
print(loo_full)
