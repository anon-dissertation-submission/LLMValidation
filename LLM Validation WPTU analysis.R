library(readxl)
library(dplyr)

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET_NO_MISSING_NO_PARTIAL <- "No Missing and No Partial"
DATA_SHEET_NO_PARTIAL <- "No Partial"
all_data <- read_excel(DATA_FILE, sheet = DATA_SHEET_NO_PARTIAL)

# Derive baseline severity per participant: earliest visit, max severity across domains
baseline_per_subj <- all_data %>%
  group_by(src_subject_id, visit) %>%
  summarise(
    visit_max_severity = max(severity, na.rm = TRUE),
    interview_date     = first(interview_date),
    .groups = "drop"
  ) %>%
  group_by(src_subject_id) %>%
  slice_min(interview_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(src_subject_id, baseline_severity = visit_max_severity)

# Determine analytic sample membership
analytic_sheet <- read_excel(DATA_FILE, sheet = DATA_SHEET_NO_MISSING_NO_PARTIAL)
analytic_ids <- unique(analytic_sheet$src_subject_id)

baseline_per_subj <- baseline_per_subj %>%
  mutate(included = src_subject_id %in% analytic_ids)

# Shapiro-Wilk
sw_included <- shapiro.test(baseline_per_subj$baseline_severity[baseline_per_subj$included])
sw_excluded <- shapiro.test(baseline_per_subj$baseline_severity[!baseline_per_subj$included])

cat("Shapiro-Wilk test for baseline severity:\n")
cat(sprintf("  Analytic sample (n=%d): W=%.3f, p=%.4f\n",
            sum(baseline_per_subj$included), sw_included$statistic, sw_included$p.value))
cat(sprintf("  Excluded sample (n=%d): W=%.3f, p=%.4f\n",
            sum(!baseline_per_subj$included), sw_excluded$statistic, sw_excluded$p.value))

# t-test and Mann-Whitney U
tt <- t.test(baseline_severity ~ included, data = baseline_per_subj)
mw <- wilcox.test(baseline_severity ~ included, data = baseline_per_subj)

cat(sprintf("\nt-test:         t=%.3f, df=%.1f, p=%.4f\n",
            tt$statistic, tt$parameter, tt$p.value))
cat(sprintf("Mann-Whitney U: W=%.1f, p=%.4f\n", mw$statistic, mw$p.value))
# ... [Your existing script, including the cat() printout statements] ...

# -------------------------------------------------------------------------
# EXPORTING TO CSV
# -------------------------------------------------------------------------

# Save the processed subject-level dataset 
# This saves the 'baseline_per_subj' dataframe with the new 'included' column
write.csv(baseline_per_subj, "baseline_severity_data.csv", row.names = FALSE)
cat("\nSaved subject data to 'baseline_severity_data.csv'\n")

# Save the statistical test results as a summary table
# First, construct a neat dataframe of your results
stats_summary <- data.frame(
  Test = c(
    "Shapiro-Wilk (Analytic Sample)", 
    "Shapiro-Wilk (Excluded Sample)", 
    "Independent t-test", 
    "Mann-Whitney U"
  ),
  Statistic = c(
    sw_included$statistic, 
    sw_excluded$statistic, 
    tt$statistic, 
    mw$statistic
  ),
  P_Value = c(
    sw_included$p.value, 
    sw_excluded$p.value, 
    tt$p.value, 
    mw$p.value
  )
)

# Then, write that summary dataframe to a CSV
write.csv(stats_summary, "baseline_severity_statistics.csv", row.names = FALSE)
cat("Saved statistical results to 'baseline_severity_statistics.csv'\n")
