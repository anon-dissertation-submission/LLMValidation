# Bland-Altman Analysis for Continuous Scoring Agreement
# LLM vs Clinician Scores

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)

# Setup & Imports

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"

df <- read_excel(DATA_FILE, sheet = DATA_SHEET)

# Data Preparation
# Ensure the required score columns are strictly numeric and drop NA rows
df_clean <- df %>%
  mutate(
    severity = as.numeric(severity),
    severity_pred = as.numeric(severity_pred),
    frequency = as.numeric(frequency),
    frequency_pred = as.numeric(frequency_pred)
  ) %>%
  drop_na(severity, severity_pred, frequency, frequency_pred)

# Bland-Altman Function
create_bland_altman <- function(data, clin_col, llm_col, metric_name) {
  
  # Extract variables
  clinician <- data[[clin_col]]
  llm <- data[[llm_col]]
  
  # Positive difference indicates systematic LLM over-scoring
  diff <- llm - clinician
  avg <- (llm + clinician) / 2
  
  # Statistical calculations
  mean_diff <- mean(diff)
  sd_diff <- sd(diff)
  
  # 95% Limits of Agreement
  loa_lower <- mean_diff - 1.96 * sd_diff
  loa_upper <- mean_diff + 1.96 * sd_diff
  
  # Percentage within ±1 point
  pct_within_1 <- mean(abs(diff) <= 1) * 100
  
  # Summary Row for CSV
  summary_row <- data.frame(
    Metric = metric_name,
    Mean_Difference = round(mean_diff, 4),
    SD_Difference = round(sd_diff, 4),
    LoA_Lower = round(loa_lower, 4),
    LoA_Upper = round(loa_upper, 4),
    Percent_Within_1_Point = round(pct_within_1, 2)
  )
  
  # Plotting dataframe
  plot_data <- data.frame(avg = avg, diff = diff)
  
  # Create plot using ggplot2 with jitter to handle discrete scoring overlap
  p <- ggplot(plot_data, aes(x = avg, y = diff)) +
    geom_jitter(width = 0.15, height = 0.15, alpha = 0.6, color = "#2166AC", size = 2, stroke = 0.5) +
    geom_hline(yintercept = mean_diff, color = "#B2182B", linetype = "solid", linewidth = 1) +
    geom_hline(yintercept = loa_upper, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = loa_lower, color = "black", linetype = "dashed", linewidth = 0.8) +
    # Add text annotations for the lines
    annotate("text", x = max(avg), y = mean_diff, label = sprintf("Mean: %.2f", mean_diff), vjust = -0.8, color = "red", hjust = 1, size = 4) +
    annotate("text", x = max(avg), y = loa_upper, label = sprintf("+1.96 SD: %.2f", loa_upper), vjust = -0.8, hjust = 1, size = 4) +
    annotate("text", x = max(avg), y = loa_lower, label = sprintf("-1.96 SD: %.2f", loa_lower), vjust = 1.5, hjust = 1, size = 4) +
    labs(
      title = paste("Bland-Altman Agreement:", metric_name, "Scores"),
      x = "Average Score (Clinician & LLM)",
      y = "Difference (LLM - Clinician)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
  
  # Save plot at 300 dpi
  file_name <- paste0("fig_bland_altman_", tolower(metric_name), ".png")
  ggsave(filename = file_name, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
  
  return(summary_row)
}

# Execute Analysis
# Process Severity and Frequency
res_severity <- create_bland_altman(df_clean, "severity", "severity_pred", "Severity")
res_frequency <- create_bland_altman(df_clean, "frequency", "frequency_pred", "Frequency")

# Formatting & Output
# Combine results and save to CSV
summary_df <- bind_rows(res_severity, res_frequency)
write_csv(summary_df, "bland_altman_summary.csv")

# Print to console to see in RStudio
print(summary_df)