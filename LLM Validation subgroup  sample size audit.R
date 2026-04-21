library(dplyr)
library(purrr)
library(readr)
library(readxl)
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
    site           = case_when(
      factor(substr(src_subject_id, 1, 2)) == "ME" ~ "ME",
      TRUE ~ "Non-ME"
    ),
    src_subject_id = as.character(src_subject_id),
    age_band       = factor(age_band),
    language       = factor(language),
    language_binary = case_when(
      tolower(trimws(language)) == "english" ~ "English",
      tolower(trimws(language)) == "unavailable" ~ "Unavailable",
      is.na(language) ~ NA_character_,
      TRUE ~ "Non-English"
    )
  )

visit_data <- collapse_to_visit_level(phenotype_data)


build_audit_rows <- function(data, demographic_label, subgroup_col) {
  data %>%
    group_by(subgroup = as.character(!!sym(subgroup_col))) %>%
    summarise(
      n_subjects   = n_distinct(src_subject_id),
      n_visits     = n(),
      n_chr_visits = sum(chr == "CHR"),
      n_hc_visits  = sum(chr == "HC"),
      .groups = "drop"
    ) %>%
    mutate(demographic = demographic_label) %>%
    select(demographic, subgroup, n_subjects, n_visits, n_chr_visits, n_hc_visits)
}

subgroup_audit <- bind_rows(
  build_audit_rows(visit_data, "Sex",       "sex"),
  build_audit_rows(visit_data, "Race",      "race"),
  build_audit_rows(visit_data, "Age band",  "age_band"),
  build_audit_rows(visit_data, "Language",  "language_binary"),
  build_audit_rows(visit_data, "Site",      "site")
)

write_csv(subgroup_audit, "subgroup_sample_sizes.csv")
print(subgroup_audit, n = Inf)
