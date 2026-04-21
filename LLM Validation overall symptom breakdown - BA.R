# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)

# Source the bootstrap helpers
source("bootstrap_helpers.R")

# Variables
N_ITERATIONS <- 10
SEED <- 13
DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"

# Data Preparation
cat("Loading and preparing data...\n")
df <- read_excel(DATA_FILE, sheet = DATA_SHEET)

df_clean <- df %>%
  drop_na(question, severity, severity_pred, frequency, frequency_pred, src_subject_id)

# Define the Balanced Accuracy Function
calc_ba <- function(truth, pred) {
  classes <- sort(unique(c(truth, pred)))
  if (length(classes) < 2) return(NA_real_)
  
  ba_vec <- numeric(length(classes))
  
  for (i in seq_along(classes)) {
    cls <- classes[i]
    
    TP <- sum(truth == cls & pred == cls, na.rm = TRUE)
    TN <- sum(truth != cls & pred != cls, na.rm = TRUE)
    FP <- sum(truth != cls & pred == cls, na.rm = TRUE)
    FN <- sum(truth == cls & pred != cls, na.rm = TRUE)
    
    sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
    spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
    
    ba_vec[i] <- (sens + spec) / 2
  }
  
  return(mean(ba_vec, na.rm = TRUE))
}

# Define the Bootstrapping Function
symptom_ba_bootstrap <- function(data, n_iterations) {
  subjects <- unique(data$src_subject_id)
  n_subj <- length(subjects)
  
  results_list <- vector("list", n_iterations)
  
  for (i in seq_len(n_iterations)) {
    sampled_ids <- sample(subjects, n_subj, replace = TRUE)
    boot_sample <- suppressMessages(get_boot_sample_subjects(data, sampled_ids))
    
    iter_res <- boot_sample %>%
      group_by(question) %>%
      summarise(
        sev_ba = calc_ba(severity, severity_pred),
        freq_ba = calc_ba(frequency, frequency_pred),
        .groups = "drop"
      ) %>%
      mutate(iteration = i)
    
    results_list[[i]] <- iter_res
    if (i %% 100 == 0) cat(sprintf("\rIteration %d/%d", i, n_iterations))
  }
  cat("\nDone!\n")
  
  summary_df <- bind_rows(results_list) %>%
    group_by(question) %>%
    summarise(
      sev_mean = mean(sev_ba, na.rm = TRUE),
      sev_lo = quantile(sev_ba, 0.025, na.rm = TRUE),
      sev_hi = quantile(sev_ba, 0.975, na.rm = TRUE),
      freq_mean = mean(freq_ba, na.rm = TRUE),
      freq_lo = quantile(freq_ba, 0.025, na.rm = TRUE),
      freq_hi = quantile(freq_ba, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_df)
}

# Run the Bootstrap
set.seed(SEED)
cat("\nRunning symptom accuracy bootstrap...\n")
boot_results <- symptom_ba_bootstrap(df_clean, N_ITERATIONS)

# Calculate exact 'n' (data points)
n_counts <- df_clean %>%
  group_by(question) %>%
  summarise(
    sev_n = sum(!is.na(severity) & !is.na(severity_pred)),
    freq_n = sum(!is.na(frequency) & !is.na(frequency_pred)),
    .groups = "drop"
  )

boot_results <- boot_results %>%
  left_join(n_counts, by = "question")

# ==============================================================================
# Format Data and Map Labels safely
# ==============================================================================
plot_data_mapped <- boot_results %>%
  mutate(
    symptom_label = case_when(
      question == "Q1" ~ "Unusual thoughts and experiences",
      question == "Q2" ~ "Suspiciousness",
      question == "Q3" ~ "Unusual somatic ideas",
      question == "Q4" ~ "Ideas of guilt",
      question == "Q5" ~ "Jealous ideas",
      question == "Q6" ~ "Unusual religious ideas",
      question == "Q7" ~ "Erotomanic ideas",
      question == "Q8" ~ "Grandiosity",
      question == "Q9" ~ "Auditory perceptual abnormalities",
      question == "Q10" ~ "Visual perceptual abnormalities",
      question == "Q11" ~ "Olfactory perceptual abnormalities",
      question == "Q12" ~ "Gustatory perceptual abnormalities",
      question == "Q13" ~ "Tactile perceptual abnormalities",
      question == "Q14" ~ "Somatic perceptual abnormalities",
      question == "Q15" ~ "Disorganised Communication Expression",
      TRUE ~ question 
    ),
    q_num = as.numeric(gsub("[^0-9]", "", question)) 
  )

ordered_labels <- plot_data_mapped %>% 
  arrange(q_num) %>% 
  pull(symptom_label) %>% 
  unique()

# ==============================================================================
# Format Data and Map Labels safely
# ==============================================================================

MIN_SUBJ <- 10

plot_data <- plot_data_mapped %>%
  pivot_longer(
    cols = c(starts_with("sev_"), starts_with("freq_")),
    names_to = c("dimension", ".value"),
    names_pattern = "(sev|freq)_(.*)"
  ) %>%
  mutate(
    dimension = ifelse(dimension == "sev", "Severity", "Frequency"),
    
    # Apply logic from Generic script for plot values and NA handling
    plot_val = ifelse(n < MIN_SUBJ, NA_real_, mean),
    
    heatmap_label = case_when(
      n < MIN_SUBJ ~ "—", 
      is.na(mean) ~ "", 
      TRUE ~ sprintf("%.2f\n(%.2f, %.2f)\nn=%d", mean, lo, hi, n)
    ),
    
    symptom_factor = factor(symptom_label, levels = rev(ordered_labels)),
    
    # Generic script text color logic (white for extremes only)
    text_color = ifelse(!is.na(plot_val) & (plot_val < 0.25 | plot_val > 0.75), "white", "black")
  ) %>% 
  ungroup()

# ==============================================================================
# Create Plot: The Heatmap
# ==============================================================================
cat("Generating heatmap...\n")

# Shared custom_theme from Generic Script
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title  = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "bottom"
  )

heatmap_plot <- ggplot(plot_data, aes(x = dimension, y = symptom_factor, fill = plot_val)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  
  # CHANGED: Increased text size from 2.8 to 3.5 to make the bold text pop!
  geom_text(aes(label = heatmap_label, colour = text_color), size = 3.5, fontface = "bold", lineheight = 0.8) +
  scale_colour_identity() + 
  
  # Exact scale configuration from Generic Script
  scale_fill_gradient2(
    low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
    midpoint = 0.50, limits = c(0, 1), na.value = "grey88", 
    name = "BA",
    guide = guide_colorbar(
      title.position = "left", title.vjust = 0.8,
      barwidth = unit(4.5, "cm"), barheight = unit(0.5, "cm"), ticks.linewidth = 0            
    )
  ) +
  labs(
    title = "Model Balanced Accuracy by Symptom (Heatmap)",
    subtitle = sprintf("Macro-Averaged Balanced Accuracy. '—' = fewer than %d subjects", MIN_SUBJ),
    x = "Rating Dimension",
    y = NULL 
  ) +
  custom_theme + 
  theme(
    plot.title.position = "plot", 
    plot.margin = margin(t = 20, r = 15, b = 15, l = 15),
    panel.grid = element_blank(),
    legend.justification = "center"
  )

fixed_width <- 7.5
fixed_height <- 11.5

# Save PDF and PNG
ggsave("symptom_ba_heatmap.pdf", plot = heatmap_plot, width = fixed_width, height = fixed_height)
ggsave("symptom_ba_heatmap.png", plot = heatmap_plot, width = fixed_width, height = fixed_height, dpi = 300)
cat(sprintf("Heatmap saved as PDF and PNG (Width: %.1f, Height: %.1f)\n", fixed_width, fixed_height))

# ==============================================================================
# Create Plot: The Point & Confidence Interval Graph
# ==============================================================================
cat("Generating CI point plot...\n")

point_plot <- ggplot(plot_data, aes(y = symptom_factor, x = mean, color = dimension, group = dimension)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), height = 0.2, linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_text(aes(x = hi + 0.02, label = sprintf("n=%d", n)), 
            position = position_dodge(width = 0.6), 
            size = 3.2, 
            show.legend = FALSE, 
            color = "gray20",
            hjust = 0) + 
  scale_color_manual(values = c("Severity" = "#08519c", "Frequency" = "#1b9e77"), name = "Dimension") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.05, 0.15))) +
  theme_minimal() +
  labs(
    title = "Model Balanced Accuracy by Symptom (95% CIs)",
    subtitle = "Macro-Averaged Balanced Accuracy. Error bars indicate 95% Confidence Intervals.",
    y = "Symptom (Question)",
    x = "Balanced Accuracy"
  ) +
  theme(
    # FIX: Align title to the far left of the entire image and add generous outer margins
    plot.title.position = "plot", 
    plot.margin = margin(t = 20, r = 15, b = 15, l = 15),
    
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 8), hjust = 0),
    plot.subtitle = element_text(size = 11, color = "gray30", margin = margin(b = 20), hjust = 0),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top"
  )

point_plot_height <- max(5.0, length(ordered_labels) * 0.8)

# ==============================================================================
# Generate and Save Tabular Summary as HTML
# ==============================================================================
cat("\nGenerating HTML summary table...\n")

# Format the data into a clean table
summary_table <- plot_data_mapped %>%
  arrange(q_num) %>%
  mutate(
    `Severity BA (95% CI)` = sprintf("%.1f%% [%.1f%%, %.1f%%]", sev_mean * 100, sev_lo * 100, sev_hi * 100),
    `Frequency BA (95% CI)` = sprintf("%.1f%% [%.1f%%, %.1f%%]", freq_mean * 100, freq_lo * 100, freq_hi * 100),
    `N (Sev/Freq)` = sprintf("%d/%d", sev_n, freq_n)
  ) %>%
  select(
    Question = question,
    Symptom = symptom_label,
    `Severity BA (95% CI)`,
    `Frequency BA (95% CI)`,
    `N (Sev/Freq)`
  )

# Create the HTML table with styling
html_table <- kable(summary_table, format = "html", escape = FALSE, align = "llccc") %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
    full_width = FALSE,
    position = "left",
    font_size = 14
  ) %>%
  add_header_above(c(" " = 2, "Balanced Accuracy Metrics" = 2, " " = 1), bold = TRUE) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#08519c") # Matches your plot theme

# Save the HTML to a file
writeLines(as.character(html_table), "symptom_ba_summary.html")
cat("Success! HTML table saved to 'symptom_ba_summary.html'\n")
