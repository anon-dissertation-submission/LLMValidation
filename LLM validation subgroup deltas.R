library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(rlang)
library(readxl)
library(tibble)
library(readr) # For write_csv

# source bootstrap_helpers.R from the same directory as this script
source("bootstrap_helpers.R")

# Settings
DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"
N_ITERATIONS <- 10   # set to 10 for a quick test, increase to 1000 for final
SEED         <- 13

# CAARMS thresholds
CAARMS_SEVERITY_LOW_APS2a   <- 3
CAARMS_SEVERITY_HIGH_APS2a  <- 5
CAARMS_FREQUENCY_LOW_APS2a  <- 3
CAARMS_FREQUENCY_HIGH_APS2a <- 6
CAARMS_SEVERITY_APS2b       <- 6
CAARMS_FREQUENCY_APS2b      <- 3
CAARMS_SEVERITY_BLIPS       <- 6
CAARMS_FREQUENCY_LOW_BLIPS  <- 4
CAARMS_FREQUENCY_HIGH_BLIPS <- 6

# Load data + derive CHR/HC labels
phenotype_data <- read_excel(DATA_FILE, sheet = DATA_SHEET) %>%
  mutate(
    sex            = factor(sex),
    race           = factor(race),
    src_subject_id = as.character(src_subject_id),
    visit          = as.character(visit), 
    age_band       = factor(age_band),
    
    # Create a unique identifier for each visit
    visit_id       = paste(src_subject_id, visit, sep = "_"),
    
    # Grouping Site: Rename ME to Melbourne
    site = factor(if_else(site == "ME", "Melbourne", "Non-Melbourne"), levels = c("Melbourne", "Non-Melbourne")),
    
    # Grouping Language into English, Unavailable, and Non-English using the 'language' column
    language_group = factor(case_when(
      language == "English" ~ "English",
      language == "Unavailable" ~ "Unavailable",
      TRUE ~ "Non-English"
    ), levels = c("English", "Non-English", "Unavailable"))
  ) %>%
  mutate(
    chr_subtype = case_when(
      severity %in% CAARMS_SEVERITY_LOW_APS2a:CAARMS_SEVERITY_HIGH_APS2a &
        frequency %in% CAARMS_FREQUENCY_LOW_APS2a:CAARMS_FREQUENCY_HIGH_APS2a ~ "APS2a",
      severity == CAARMS_SEVERITY_APS2b &
        frequency == CAARMS_FREQUENCY_APS2b ~ "APS2b",
      severity == CAARMS_SEVERITY_BLIPS &
        frequency %in% CAARMS_FREQUENCY_LOW_BLIPS:CAARMS_FREQUENCY_HIGH_BLIPS ~ "BLIPS",
      TRUE ~ "Below threshold"
    ),
    chr_subtype_pred = case_when(
      severity_pred %in% CAARMS_SEVERITY_LOW_APS2a:CAARMS_SEVERITY_HIGH_APS2a &
        frequency_pred %in% CAARMS_FREQUENCY_LOW_APS2a:CAARMS_FREQUENCY_HIGH_APS2a ~ "APS2a",
      severity_pred == CAARMS_SEVERITY_APS2b &
        frequency_pred == CAARMS_FREQUENCY_APS2b ~ "APS2b",
      severity_pred == CAARMS_SEVERITY_BLIPS &
        frequency_pred %in% CAARMS_FREQUENCY_LOW_BLIPS:CAARMS_FREQUENCY_HIGH_BLIPS ~ "BLIPS",
      TRUE ~ "Below threshold"
    ),
    chr      = factor(if_else(chr_subtype      == "Below threshold", "HC", "CHR"), levels = c("HC", "CHR")),
    chr_pred = factor(if_else(chr_subtype_pred == "Below threshold", "HC", "CHR"), levels = c("HC", "CHR"))
  )

vars_needed <- c("src_subject_id", "visit", "visit_id", "chr", "chr_pred", "sex", "age_band", "race", "site", "language_group")
# Global filter for complete cases
df <- phenotype_data %>% drop_na(all_of(vars_needed))

cat(sprintf("Rows after drop_na: %d | Unique subjects: %d | Unique visits: %d\n", 
            nrow(df), n_distinct(df$src_subject_id), n_distinct(df$visit_id)))

# Cluster bootstrap helper
cluster_boot_sample <- function(data, vis_ids) {
  sampled <- sample(vis_ids, size = length(vis_ids), replace = TRUE)
  tbl     <- tibble(visit_id = sampled, .rep = seq_along(sampled))
  data %>%
    inner_join(tbl, by = "visit_id", relationship = "many-to-many") %>%
    select(-.rep)
}

# Bootstrap ALL metrics for one grouping variable
boot_fairness_onevar <- function(data, group_var, ref_group,
                                 n_boot   = N_ITERATIONS,
                                 seed     = SEED,
                                 min_vis  = 10) { 
  
  set.seed(seed)
  gsym <- sym(group_var)
  
  d0 <- data %>% filter(!is.na(!!gsym)) %>% mutate(.g = as.character(!!gsym))
  available_groups <- unique(d0$.g)
  
  if (!ref_group %in% available_groups) {
    stop(sprintf("Reference group '%s' not found. Available: %s", ref_group, paste(sort(available_groups), collapse = ", ")))
  }
  
  denom <- d0 %>% group_by(.g) %>% summarise(n_visits = n_distinct(visit_id), .groups = "drop")
  keep_groups <- denom %>% filter(n_visits >= min_vis) %>% pull(.g)
  
  if (!ref_group %in% keep_groups) {
    stop(sprintf("Reference group '%s' has fewer than %d visits.", ref_group, min_vis))
  }
  
  # EXPLICITLY CAPTURE THE REFERENCE GROUP'S N
  ref_n <- denom %>% filter(.g == ref_group) %>% pull(n_visits)
  
  d0 <- d0 %>% filter(.g %in% keep_groups)
  comp_groups <- sort(setdiff(unique(d0$.g), ref_group))
  
  if (length(comp_groups) == 0) return(tibble())
  
  all_vis_ids <- unique(d0$visit_id)
  
  cat(sprintf("\n[%s] ref='%s' | compare: %s | visits: %d | boots: %d\n",
              group_var, ref_group, paste(comp_groups, collapse = ", "), length(all_vis_ids), n_boot))
  
  # Storage arrays
  boot_ba   <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  boot_tpr  <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  boot_fpr  <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  
  for (b in seq_len(n_boot)) {
    boot_df <- cluster_boot_sample(d0, all_vis_ids)
    
    # Calculate BA, TPR (sens), and Spec for every group
    metrics_by_g <- boot_df %>%
      group_by(.g) %>%
      group_modify(~ {
        m <- calc_metrics(.x)
        tibble(bal_acc = m$bal_acc, sens = m$sens, spec = m$spec)
      }) %>%
      ungroup()
    
    ref_row <- metrics_by_g %>% filter(.g == ref_group)
    if (nrow(ref_row) != 1 || is.na(ref_row$bal_acc)) next
    
    for (g in comp_groups) {
      g_row <- metrics_by_g %>% filter(.g == g)
      if (nrow(g_row) == 1 && !is.na(g_row$bal_acc)) {
        boot_ba[b, g]  <- g_row$bal_acc - ref_row$bal_acc
        boot_tpr[b, g] <- g_row$sens - ref_row$sens
        # delta FPR = -(delta Specificity)
        boot_fpr[b, g] <- -(g_row$spec - ref_row$spec) 
      }
    }
    if (b %% 100 == 0) cat(sprintf("  iter %d/%d\n", b, n_boot))
  }
  
  # Summarize for each metric
  summarise_metric <- function(mat, metric_name) {
    lapply(comp_groups, function(g) {
      vals <- mat[, g]
      vals <- vals[!is.na(vals)]
      
      # Calculate empirical p-value (two-sided)
      p_empirical <- 2 * min(mean(vals >= 0), mean(vals <= 0))
      
      tibble(
        .g = g,
        metric = metric_name,
        estimate = mean(vals),
        lo = quantile(vals, 0.025),
        hi = quantile(vals, 0.975),
        p_val = p_empirical,
        n_valid = length(vals)
      )
    }) %>% bind_rows()
  }
  
  bind_rows(
    summarise_metric(boot_ba, "\u0394 Balanced Accuracy"),
    summarise_metric(boot_tpr, "\u0394 TPR (Sensitivity)"),
    summarise_metric(boot_fpr, "\u0394 FPR (1 - Specificity)")
  ) %>%
    left_join(denom, by = ".g") %>%
    mutate(group_var = group_var, ref = ref_group, ref_n_visits = ref_n)
}

# Run bootstrap for all subgroups
d_sex  <- boot_fairness_onevar(df, "sex",            ref_group = "F")
d_age  <- boot_fairness_onevar(df, "age_band",       ref_group = "18-25")
d_race <- boot_fairness_onevar(df, "race",           ref_group = "White")
d_site <- boot_fairness_onevar(df, "site",           ref_group = "Melbourne") 
d_lang <- boot_fairness_onevar(df, "language_group", ref_group = "English") 

# Combine and prepare plotting table
plot_df <- bind_rows(
  d_sex  %>% mutate(panel = "Sex"),
  d_age  %>% mutate(panel = "Age band"),
  d_race %>% mutate(panel = "Race"),
  d_site %>% mutate(panel = "Site"),
  d_lang %>% mutate(panel = "Language")
) %>%
  filter(!is.na(estimate)) %>%
  mutate(
    # Compute FDR (Benjamini-Hochberg) across all tests
    fdr = p.adjust(p_val, method = "BH"),
    # Create a flag for significance
    fdr_sig = factor(if_else(fdr < 0.05, "Significant (FDR < 0.05)", "Not Significant"), 
                     levels = c("Significant (FDR < 0.05)", "Not Significant")),
    
    panel = factor(panel, levels = c("Sex", "Age band", "Race", "Site", "Language")),
    metric = factor(metric, levels = c("\u0394 Balanced Accuracy", "\u0394 TPR (Sensitivity)", "\u0394 FPR (1 - Specificity)")),
    label = paste0(.g, " vs ", ref, "\n(n=", n_visits, ", ref_n=", ref_n_visits, ")")
  ) %>%
  group_by(panel, metric) %>%
  mutate(label = fct_reorder(label, estimate)) %>%
  ungroup()

# Export to CSV
csv_out <- plot_df %>%
  select(panel, group_var, 
         ref_group = ref, reference_n_visits = ref_n_visits, 
         comparison_group = .g, comparison_n_visits = n_visits, 
         metric, estimate, lower_ci_95 = lo, upper_ci_95 = hi, 
         empirical_p_val = p_val, fdr) %>%
  arrange(panel, metric, estimate)

write_csv(csv_out, "bootstrap_fairness_metrics_subgroups_by_visit.csv")
cat("\nSaved data to: bootstrap_fairness_metrics_subgroups_by_visit.csv\n")

# Multi-Panel Forest Plot
# Removed alpha mapping so all lines and points are fully opaque
p_forest <- ggplot(plot_df, aes(x = estimate, y = label, colour = panel, shape = fdr_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  # Error bars are now fully opaque regardless of significance
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.8) +
  # Plot points map shape to significance (Solid vs Open)
  geom_point(size = 3) +
  scale_colour_brewer(palette = "Dark2") +
  # Map Shapes: 16 = Solid Circle, 1 = Open Circle
  scale_shape_manual(values = c("Significant (FDR < 0.05)" = 16, "Not Significant" = 1), name = "Significance") +
  facet_grid(panel ~ metric, scales = "free_y", space = "free_y") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Visit-level binary fairness across demographic subgroups", 
    x        = "Difference relative to reference group (95% cluster-bootstrap CI)",
    y        = NULL,
    caption  = sprintf(
      "Cluster bootstrap by visit (n = %d iterations).\n\u0394 FPR calculated as -( \u0394 Specificity ).\nSignificance adjusted using Benjamini-Hochberg FDR.",
      N_ITERATIONS
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text.x       = element_text(face = "bold", size = 11),
    strip.text.y       = element_text(face = "bold", size = 11, angle = 270),
    axis.text.y        = element_text(face = "bold", size = 10, lineheight = 1.1), 
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.caption       = element_text(size = 9, colour = "grey50"),
    # Keep the significance legend at the bottom, but hide the panel color legend
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold")
  ) +
  guides(colour = "none") # Hides the redundant color legend

print(p_forest)

ggsave("forest_fairness_metrics_grid_by_visit.png", p_forest, width = 12, height = 11, dpi = 300)
cat("\nSaved plot to: forest_fairness_metrics_grid_by_visit.png\n")