# Continuous Scoring Assessment: LLM vs Clinician Scores (Parallelized)

library(dplyr)
library(tidyr)
library(readr)
library(irr)
library(tibble)
library(readxl)
library(purrr)
library(future)
library(furrr)
library(progressr) # Added for parallel progress tracking

# Setup & Imports
# Adjust the data file path if your data is located elsewhere
DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"
N_BOOT     <- 10000  # Standard for stable CIs
SEED       <- 13

set.seed(SEED)

# Load data
df <- read_excel(DATA_FILE, sheet = DATA_SHEET)

# Create a unique visit ID to bootstrap at the visit level
df <- df %>% mutate(visit_id = paste(src_subject_id, visit, sep = "_"))

# Ensure the required score columns are strictly numeric
df <- df %>%
  mutate(
    severity = as.numeric(severity),
    severity_pred = as.numeric(severity_pred),
    frequency = as.numeric(frequency),
    frequency_pred = as.numeric(frequency_pred)
  )

# Metrics Calculation Function
calc_continuous_metrics <- function(data, true_col, pred_col) {
  d_clean <- data %>%
    select(all_of(c(true_col, pred_col))) %>%
    drop_na()
  
  if(nrow(d_clean) < 2) return(c(icc = NA_real_, r = NA_real_, mae = NA_real_))
  
  r_val <- cor(d_clean[[true_col]], d_clean[[pred_col]], method = "pearson")
  mae_val <- mean(abs(d_clean[[true_col]] - d_clean[[pred_col]]))
  
  icc_res <- irr::icc(d_clean, model = "twoway", type = "agreement", unit = "single")
  icc_val <- icc_res$value
  
  return(c(icc = icc_val, r = r_val, mae = mae_val))
}

# Point Estimates
cat("Calculating point estimates...\n")
obs_sev  <- calc_continuous_metrics(df, "severity", "severity_pred")
obs_freq <- calc_continuous_metrics(df, "frequency", "frequency_pred")

# Visit-Level Bootstrapping (PARALLEL with PROGRESS)
cat(sprintf("Running %d bootstrap iterations at the visit level in PARALLEL...\n", N_BOOT))

unique_visits <- unique(df$visit_id)

# Set up parallel backend
plan(multisession, workers = availableCores())

# Enable progress handlers globally in your R session
handlers(global = TRUE)
handlers("progress") # Uses the standard progress bar format

# Wrap the execution in with_progress()
with_progress({
  
  # Initialize the progressor
  p <- progressor(steps = N_BOOT)
  
  # Parallel Bootstrap Loop
  boot_res <- future_map_dfr(seq_len(N_BOOT), function(i) {
    
    # 1. Sample visits with replacement
    sampled_visits <- sample(unique_visits, replace = TRUE)
    
    # 2. Reconstruct dataset efficiently using a frequency table join
    boot_counts <- as.data.frame(table(visit_id = sampled_visits), stringsAsFactors = FALSE)
    boot_df <- df %>%
      inner_join(boot_counts, by = "visit_id", relationship = "many-to-many") %>%
      uncount(Freq)
    
    # 3. Calculate metrics
    res_sev <- calc_continuous_metrics(boot_df, "severity", "severity_pred")
    res_freq <- calc_continuous_metrics(boot_df, "frequency", "frequency_pred")
    
    # 4. Signal progress back to the main session
    p()
    
    # 5. Return as a single-row tibble
    tibble(
      icc_sev = res_sev["icc"],
      r_sev   = res_sev["r"],
      mae_sev = res_sev["mae"],
      icc_freq = res_freq["icc"],
      r_freq   = res_freq["r"],
      mae_freq = res_freq["mae"]
    )
    
  }, .options = furrr_options(seed = TRUE, packages = c("dplyr", "tidyr", "irr")))
  
})

# Formatting & Output
format_stat <- function(obs, boot_vec) {
  v <- na.omit(boot_vec)
  if (length(v) == 0 || is.na(obs)) return("NA")
  
  lower <- quantile(v, 0.025)
  upper <- quantile(v, 0.975)
  sprintf("%.3f [%.3f, %.3f]", obs, lower, upper)
}

final_table <- tibble(
  Metric = c(
    "Intraclass Correlation (ICC) [95% CI]",
    "Pearson's r [95% CI]",
    "Mean Absolute Error (MAE) [95% CI]"
  ),
  `Severity Score [95% CI]` = c(
    format_stat(obs_sev["icc"], boot_res$icc_sev),
    format_stat(obs_sev["r"], boot_res$r_sev),
    format_stat(obs_sev["mae"], boot_res$mae_sev)
  ),
  `Frequency Score [95% CI]` = c(
    format_stat(obs_freq["icc"], boot_res$icc_freq),
    format_stat(obs_freq["r"], boot_res$r_freq),
    format_stat(obs_freq["mae"], boot_res$mae_freq)
  )
)

cat("\n=======================================================================\n")
cat("Required output — Overall Continuous Scoring\n")
cat("=======================================================================\n\n")

print(knitr::kable(final_table, align = "lcc"))

write_csv(final_table, "continuous_scoring_metrics.csv")
