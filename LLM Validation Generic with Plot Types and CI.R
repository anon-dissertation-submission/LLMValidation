# Generic validation script with flexible plot types and confidence interval generation.

library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(rlang)
library(boot)
library(irr)      
library(Metrics)  

N_BOOT   <- 10
SEED     <- 13   
MIN_SUBJ <- 10    

predictor_categories <- list(
  sex      = "Female (F)",
  age_band = "18-24",            
  language = "English"  
)


# DATA PREPARATION 

CAARMS_SEVERITY_LOW_APS2a <- 3
CAARMS_SEVERITY_HIGH_APS2a <- 5
CAARMS_FREQUENCY_LOW_APS2a <- 3
CAARMS_FREQUENCY_HIGH_APS2a <- 6
CAARMS_SEVERITY_APS2b <- 6
CAARMS_FREQUENCY_APS2b <- 3
CAARMS_SEVERITY_BLIPS <- 6
CAARMS_FREQUENCY_LOW_BLIPS <- 4
CAARMS_FREQUENCY_HIGH_BLIPS <- 6

DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"

# CAARMS criteria for determining if a symptom is in the CHR domain (binary classification)
is_chr_domain <- function(sev, freq) {
  (sev %in% CAARMS_SEVERITY_LOW_APS2a:CAARMS_SEVERITY_HIGH_APS2a &
     freq %in% CAARMS_FREQUENCY_LOW_APS2a:CAARMS_FREQUENCY_HIGH_APS2a) |
    (sev == CAARMS_SEVERITY_APS2b  & freq == CAARMS_FREQUENCY_APS2b) |
    (sev == CAARMS_SEVERITY_BLIPS  & freq %in% CAARMS_FREQUENCY_LOW_BLIPS:CAARMS_FREQUENCY_HIGH_BLIPS)
}

# Domain mapping
domain_map <- tibble::tribble(
  ~q_num, ~domain_lab,
  1L, "Unusual thoughts and experiences",
  2L, "Suspiciousness",
  3L, "Unusual somatic ideas",
  4L, "Ideas of guilt",
  5L, "Jealous ideas",
  6L, "Unusual religious ideas",
  7L, "Erotomanic ideas",
  8L, "Grandiosity",
  9L, "Auditory perceptual abnormalities",
  10L, "Visual perceptual abnormalities",
  11L, "Olfactory perceptual abnormalities",
  12L, "Gustatory perceptual abnormalities",
  13L, "Tactile perceptual abnormalities",
  14L, "Somatic perceptual abnormalities",
  15L, "Disorganised Communication Expression"
)

# Read in the raw data and prepare it for analysis
raw_data <- read_excel(DATA_FILE, sheet = DATA_SHEET)

raw <- raw_data %>%
  mutate(
    src_subject_id = as.character(src_subject_id),
    sex      = factor(sex),
    age_band = factor(age_band)
  )

if ("language" %in% colnames(raw)) {
  raw <- raw %>%
    mutate(
      language = case_when(
        is.na(language) ~ NA_character_,
        tolower(trimws(language)) == "english" ~ "English",
        TRUE ~ "Non-English"
      )
    )
} 

domain_rows <- raw %>%
  mutate(q_num = readr::parse_number(as.character(question))) %>%
  left_join(domain_map, by = "q_num") %>%
  filter(
    !is.na(severity),      !is.na(frequency),
    !is.na(severity_pred), !is.na(frequency_pred)
  ) %>%
  mutate(
    chr_true = as.integer(is_chr_domain(severity,      frequency)),
    chr_pred = as.integer(is_chr_domain(severity_pred, frequency_pred))
  )


# Metrics and Bootstrapping Helper functions
calc_class_metrics <- function(truth, pred) {
  truth <- as.integer(truth); pred  <- as.integer(pred)
  n_pos <- sum(truth == 1, na.rm = TRUE); n_neg <- sum(truth == 0, na.rm = TRUE)
  
  tpr <- if (n_pos > 0) sum(pred == 1 & truth == 1, na.rm = TRUE) / n_pos else NA_real_
  tnr <- if (n_neg > 0) sum(pred == 0 & truth == 0, na.rm = TRUE) / n_neg else NA_real_
  fpr <- if (n_neg > 0) sum(pred == 1 & truth == 0, na.rm = TRUE) / n_neg else NA_real_
  ba  <- mean(c(tpr, tnr), na.rm = TRUE)
  
  c(BA = ba, TPR = tpr, TNR = tnr, FPR = fpr)
}

calc_cont_metrics <- function(truth_score, pred_score) {
  valid <- !is.na(truth_score) & !is.na(pred_score)
  ts <- truth_score[valid]; ps <- pred_score[valid]
  if (length(ts) < 2) return(c(ICC = NA_real_, Pearson = NA_real_, MAE = NA_real_))
  
  r <- cor(ts, ps, method = "pearson")
  mae_val <- mae(ts, ps)
  icc_res <- tryCatch({
    icc(data.frame(ts, ps), model = "twoway", type = "agreement", unit = "single")$value
  }, error = function(e) NA_real_)
  
  c(ICC = icc_res, Pearson = r, MAE = mae_val)
}

# Main cluster bootstrap function to sample subjects and compute metrics by group and domain
cluster_boot_sample <- function(data, subj_ids) {
  sampled <- sample(subj_ids, size = length(subj_ids), replace = TRUE)
  tbl     <- tibble(src_subject_id = sampled, .rep = seq_along(sampled))
  data %>% inner_join(tbl, by = "src_subject_id", relationship = "many-to-many") %>% select(-.rep)
}

boot_all_metrics_by_domain <- function(data, group_var, n_boot = N_BOOT, seed = SEED) {
  set.seed(seed)
  gsym <- rlang::sym(group_var)
  
  d0 <- data %>% filter(!is.na(!!gsym)) %>% mutate(.grp = as.character(!!gsym))
  groups  <- sort(unique(d0$.grp)); domains <- sort(unique(d0$q_num))
  
  subj_per_cell <- d0 %>% group_by(.grp, q_num) %>% summarise(n_subj = n_distinct(src_subject_id), .groups = "drop")
  all_ids <- unique(d0$src_subject_id)
  
  metric_names <- c("BA", "TPR", "TNR", "FPR", "ICC", "Pearson", "MAE")
  boot_store <- list()
  for (m in metric_names) {
    boot_store[[m]] <- lapply(groups, function(g) matrix(NA_real_, nrow = n_boot, ncol = length(domains), dimnames = list(NULL, as.character(domains))))
    names(boot_store[[m]]) <- groups
  }
  
  cat(sprintf("\n[%s] subjects: %d | boots: %d\n", group_var, length(all_ids), n_boot))
  
  for (b in seq_len(n_boot)) {
    boot_df <- cluster_boot_sample(d0, all_ids)
    for (g in groups) {
      g_df <- boot_df %>% filter(.grp == g)
      if (nrow(g_df) == 0) next
      for (qn in domains) {
        q_df <- g_df %>% filter(q_num == qn)
        if (nrow(q_df) < 2) next
        c_mets <- calc_class_metrics(q_df$chr_true, q_df$chr_pred)
        co_mets <- calc_cont_metrics(q_df$severity, q_df$severity_pred)
        
        boot_store[["BA"]][[g]][b, as.character(qn)]      <- c_mets["BA"]
        boot_store[["TPR"]][[g]][b, as.character(qn)]     <- c_mets["TPR"]
        boot_store[["TNR"]][[g]][b, as.character(qn)]     <- c_mets["TNR"]
        boot_store[["FPR"]][[g]][b, as.character(qn)]     <- c_mets["FPR"]
        boot_store[["ICC"]][[g]][b, as.character(qn)]     <- co_mets["ICC"]
        boot_store[["Pearson"]][[g]][b, as.character(qn)] <- co_mets["Pearson"]
        boot_store[["MAE"]][[g]][b, as.character(qn)]     <- co_mets["MAE"]
      }
    }
    if (b %% 100 == 0) cat(sprintf("  iter %d / %d\n", b, n_boot))
  }
  
  final_res <- lapply(metric_names, function(m) {
    lapply(groups, function(g) {
      mat <- boot_store[[m]][[g]]
      lapply(as.character(domains), function(qn) {
        vals <- mat[, qn]; vals_ok <- vals[!is.na(vals)]
        n_s <- subj_per_cell %>% filter(.grp == g, q_num == as.integer(qn)) %>% pull(n_subj)
        if (length(n_s) == 0) n_s <- 0
        tibble(
          Metric = m, group_var = group_var, group = g, q_num = as.integer(qn), n_subj = n_s,
          mean_val = if(length(vals_ok) > 0) mean(vals_ok) else NA_real_,
          ci_lo    = if(length(vals_ok) > 1) quantile(vals_ok, 0.025) else NA_real_,
          ci_hi    = if(length(vals_ok) > 1) quantile(vals_ok, 0.975) else NA_real_
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  final_res %>% left_join(domain_map, by = "q_num")
}


# Plotting functions for Heatmap, Dot-and-Whisker, and Bar Chart with CI and N labels
# Shared theme
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title  = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "bottom"
  )

# Updated Heatmap: Shows Mean, CI [lo, hi], and N in every tile
plot_heatmap <- function(df_summary, metric_filter) {
  # Note: 'Metric' must be capitalized to match your bootstrap output
  plot_data <- df_summary %>% filter(Metric == metric_filter) %>%
    mutate(
      plot_val = ifelse(n_subj < MIN_SUBJ, NA_real_, mean_val),
      # This creates the multi-line label for the tile
      tile_label = case_when(
        n_subj < MIN_SUBJ ~ "—", 
        is.na(mean_val) ~ "", 
        TRUE ~ sprintf("%.2f\n[%.2f, %.2f]\nn=%d", mean_val, ci_lo, ci_hi, n_subj)
      ),
      txt_col = ifelse(!is.na(plot_val) & (plot_val < 0.25 | plot_val > 0.75), "white", "black")
    )
  
  ggplot(plot_data, aes(x = group, y = reorder(domain_lab, desc(q_num)), fill = plot_val)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = tile_label, colour = txt_col), size = 2.8, fontface = "bold", lineheight = 0.8) +
    scale_colour_identity() +
    scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
                         midpoint = 0.50, limits = c(0, 1), na.value = "grey88", 
                         name = paste(metric_filter)) +
    labs(
      title = sprintf("%s by Symptom Domain and %s", metric_filter, unique(plot_data$group_var)),
      subtitle = sprintf("Tiles show Mean [95%% CI] and Sample Size (n). '—' = < %d subjects", MIN_SUBJ), 
      x = NULL, y = NULL
    ) + custom_theme + theme(panel.grid = element_blank())
}

# Updated Dot-and-Whisker: Adds text label for N next to each point
plot_grouped_point <- function(df_summary, metric_filter) {
  plot_data <- df_summary %>% filter(Metric == metric_filter, n_subj >= MIN_SUBJ)
  
  ggplot(plot_data, aes(x = reorder(domain_lab, desc(q_num)), y = mean_val, color = group)) +
    # Error bars for the 95% CI
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), 
                    position = position_dodge(width = 0.6), size = 0.6) +
    # Text label for sample size (n)
    geom_text(aes(label = paste0("n=", n_subj)), 
              position = position_dodge(width = 0.6), vjust = -1.2, size = 3, show.legend = FALSE) +
    coord_flip() +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = sprintf("%s by Symptom Domain and %s", metric_filter, unique(plot_data$group_var)),
      subtitle = "Whiskers show 95% Confidence Interval; n indicates sample size", 
      x = NULL, y = metric_filter, color = "Subgroup:"
    ) + custom_theme + theme(panel.grid.major.y = element_blank())
}

# Updated Bar Chart: Adds text label for N on top of each bar
plot_grouped_bar <- function(df_summary, metric_filter) {
  plot_data <- df_summary %>% filter(Metric == metric_filter, n_subj >= MIN_SUBJ)
  
  ggplot(plot_data, aes(x = reorder(domain_lab, desc(q_num)), y = mean_val, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    # Error bars for the 95% CI
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), 
                  position = position_dodge(width = 0.7), width = 0.25, color = "black") +
    # Text label for sample size (n)
    geom_text(aes(label = paste0("n=", n_subj)), 
              position = position_dodge(width = 0.7), vjust = -0.5, size = 3, show.legend = FALSE) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = sprintf("%s by Symptom Domain and %s", metric_filter, unique(plot_data$group_var)),
      subtitle = "Error bars show 95% Confidence Interval; n indicates sample size", 
      x = NULL, y = metric_filter, fill = "Subgroup:"
    ) + custom_theme + theme(panel.grid.major.y = element_blank())
}

# Ovarall accuracy heatmap and bar chart (no subgroups, just overall metrics)

# Create a dummy "Overall" category
domain_rows_overall <- domain_rows %>% 
  mutate(Overall = "All Patients")

# Run the bootstrap analysis for the Overall category
summary_overall <- boot_all_metrics_by_domain(domain_rows_overall, "Overall")

# Generate the General Accuracy Heatmap (Balanced Accuracy)
p_ba_overall <- plot_heatmap(summary_overall, "BA") +
  labs(
    title = "General Model Accuracy (BA) by Symptom Domain",
    subtitle = "Analysis across all participants (No subgroups)",
    x = NULL
  )

# Save the plot
ggsave("plot_BA_Overall_heatmap.png", p_ba_overall, width = 6, height = 8, dpi = 300)

# (Optional) Generate a Bar Chart for General Accuracy
# This is often easier to read for a single category
p_ba_bar <- summary_overall %>%
  filter(Metric == "BA", n_subj >= MIN_SUBJ) %>%
  ggplot(aes(x = reorder(domain_lab, mean_val), y = mean_val)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_text(aes(label = sprintf("%.2f (n=%d)", mean_val, n_subj)), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  ylim(0, 1.1) +
  labs(title = "General Model Accuracy by Symptom", 
       subtitle = "Bars show Mean Balanced Accuracy with 95% CI",
       x = "Symptom Domain", y = "Balanced Accuracy") +
  custom_theme

ggsave("plot_BA_Overall_bar.png", p_ba_bar, width = 8, height = 7, dpi = 300)


