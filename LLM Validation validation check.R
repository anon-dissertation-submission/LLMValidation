source("bootstrap_helpers.R")
source("visit_level_helpers.R")
library(readxl)

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"
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

# Validation checks
stopifnot(nrow(visit_data) == 547)
stopifnot(length(unique(visit_data$src_subject_id)) == 331)

metrics <- calc_metrics(visit_data)
print(metrics)
# Expected: sens ≈ 0.933, spec ≈ 0.576, bal_acc ≈ 0.755
write.csv(metrics, "validation_check.csv", row.names = FALSE)
cat("Saved statistical results to 'validation_check.csv'\n")


boundary_check <- visit_data %>%
  filter(race %in% c("American Indian/Alaska Native",
                     "More than one race",
                     "Unknown or not reported",
                     "White")) %>%
  group_by(race) %>%
  summarise(
    n_visits              = n(),
    n_subjects            = n_distinct(src_subject_id),
    n_chr_positive_visits = sum(chr == "CHR"),
    n_true_positives      = sum(chr == "CHR" & chr_pred == "CHR"),
    n_false_negatives     = sum(chr == "CHR" & chr_pred == "HC"),
    TPR                   = n_true_positives / n_chr_positive_visits,
    .groups = "drop"
  ) %>%
  arrange(match(race, c("White", "American Indian/Alaska Native",
                        "More than one race", "Unknown or not reported")))

write.csv(boundary_check, "boundary_artefact_check.csv")
cat("Saved boundary artefact check to 'boundary_artefact_check.csv'\n")
print(boundary_check)
