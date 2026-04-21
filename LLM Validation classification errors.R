# Error Analysis: Domain-Level Drivers of Visit-Level Misclassifications

library(dplyr)
library(tidyr)
library(readxl)
library(readr)

# Setup & Imports
source("bootstrap_helpers.R")

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"

# Domain Mapping
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

# Classification Logic
is_chr_domain <- function(sev, freq) {
  (sev %in% 3:5 & freq %in% 3:6) | # APS2a
    (sev == 6 & freq == 3) |         # APS2b
    (sev == 6 & freq %in% 4:6)       # BLIPS
}

# Load Data and Prepare Observations
cat("Loading data and classifying observations...\n")

raw_data <- read_excel(DATA_FILE, sheet = DATA_SHEET) %>%
  rename(q_num = question) %>%
  mutate(src_subject_id = as.character(src_subject_id),
         visit = as.character(visit))

obs_df <- raw_data %>%
  filter(!is.na(severity), !is.na(frequency), !is.na(severity_pred), !is.na(frequency_pred)) %>%
  left_join(domain_map, by = "q_num") %>%
  filter(!is.na(domain_lab)) %>%
  mutate(
    chr_true_bin = if_else(is_chr_domain(severity, frequency), 1, 0),
    chr_pred_bin = if_else(is_chr_domain(severity_pred, frequency_pred), 1, 0)
  )

# Aggregate to Visit Level and Identify Misclassifications
cat("Aggregating to visit level...\n")

visit_df <- obs_df %>%
  group_by(src_subject_id, visit) %>%
  summarise(
    visit_true_chr = if_else(any(chr_true_bin == 1), 1, 0),
    visit_pred_chr = if_else(any(chr_pred_bin == 1), 1, 0),
    .groups = "drop"
  ) %>%
  mutate(
    is_fp = if_else(visit_pred_chr == 1 & visit_true_chr == 0, 1, 0),
    is_fn = if_else(visit_pred_chr == 0 & visit_true_chr == 1, 1, 0)
  )

# Total counts of misclassified visits
total_fp_visits <- sum(visit_df$is_fp)
total_fn_visits <- sum(visit_df$is_fn)

cat(sprintf("Found %d False Positive visits and %d False Negative visits.\n", 
            total_fp_visits, total_fn_visits))

# Analyze False Positives
# Which domains had chr_pred == 1 in visits that were actually HC?
cat("Analyzing False Positive drivers...\n")

fp_visit_ids <- visit_df %>% filter(is_fp == 1) %>% select(src_subject_id, visit)

fp_drivers <- obs_df %>%
  inner_join(fp_visit_ids, by = c("src_subject_id", "visit")) %>%
  group_by(Domain = domain_lab) %>%
  summarise(
    n_fp_visits = sum(chr_pred_bin == 1),
    .groups = "drop"
  ) %>%
  mutate(
    `Total FP visits` = total_fp_visits,
    `% FP visits where domain drove error` = (n_fp_visits / total_fp_visits) * 100
  ) %>%
  arrange(desc(n_fp_visits)) %>%
  mutate(Rank = row_number())

write_csv(fp_drivers, "domain_cooccurrence_false_positives.csv")

# Analyze False Negatives
# Which domains had chr_true == 1 but model missed them (chr_pred == 0)?
cat("Analyzing False Negative drivers...\n")

fn_visit_ids <- visit_df %>% filter(is_fn == 1) %>% select(src_subject_id, visit)

fn_drivers <- obs_df %>%
  inner_join(fn_visit_ids, by = c("src_subject_id", "visit")) %>%
  group_by(Domain = domain_lab) %>%
  summarise(
    n_fn_visits = sum(chr_true_bin == 1 & chr_pred_bin == 0),
    .groups = "drop"
  ) %>%
  mutate(
    `Total FN visits` = total_fn_visits,
    `% FN visits where domain was missed` = (n_fn_visits / total_fn_visits) * 100
  ) %>%
  arrange(desc(n_fn_visits)) %>%
  mutate(Rank = row_number())

write_csv(fn_drivers, "domain_cooccurrence_false_negatives.csv")

# Final Display
cat("\n--- Top Drivers of False Positives ---\n")
print(head(fp_drivers, 5))

cat("\n--- Top Domains Missed in False Negatives ---\n")
print(head(fn_drivers, 5))

cat("\nAnalysis complete. CSV files saved.\n")