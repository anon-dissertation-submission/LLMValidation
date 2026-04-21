# LLM Validation: Continuous Metrics Subgroup Analysis with Cluster Bootstrap (Parallelized)
# This script performs a cluster bootstrap analysis to evaluate the differences in continuous performance metrics
# (ICC, MAE, Pearson's r) across demographic subgroups. The bootstrap is parallelized for efficiency.

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(rlang)
library(readxl)
library(tibble)
library(readr)
library(irr)
library(future)
library(future.apply)
library(progressr)

source("bootstrap_helpers.R")

# Settings
DATA_FILE  <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET <- "No Missing and No Partial"
N_ITERATIONS <- 10   # set to 10 for a quick test, increase to 10000 for final
SEED         <- 13

# Configure the progress bar to show in the console
handlers(global = TRUE)
handlers("txtprogressbar")

# Continuous Metrics Calculator
calc_continuous_metrics <- function(true_vals, pred_vals) {
  valid <- !is.na(true_vals) & !is.na(pred_vals)
  t <- true_vals[valid]
  p <- pred_vals[valid]
  
  if (length(t) < 3) {
    return(list(icc = NA_real_, mae = NA_real_, pearson = NA_real_))
  }
  
  mae <- mean(abs(t - p))
  pearson <- suppressWarnings(cor(t, p, method = "pearson"))
  
  icc_val <- tryCatch({
    res <- irr::icc(cbind(t, p), model = "twoway", type = "agreement", unit = "single")
    res$value
  }, error = function(e) NA_real_)
  
  list(icc = icc_val, mae = mae, pearson = pearson)
}

# Load data + derive continuous variables
phenotype_data <- read_excel(DATA_FILE, sheet = DATA_SHEET) %>%
  mutate(
    sex            = factor(sex),
    race           = factor(race),
    src_subject_id = as.character(src_subject_id),
    visit          = as.character(visit), 
    age_band       = factor(age_band),
    severity       = as.numeric(severity),
    severity_pred  = as.numeric(severity_pred),
    frequency      = as.numeric(frequency),
    frequency_pred = as.numeric(frequency_pred),
    visit_id       = paste(src_subject_id, visit, sep = "_"),
    site = factor(if_else(site == "ME", "Melbourne", "Non-Melbourne"), levels = c("Melbourne", "Non-Melbourne")),
    language_group = factor(case_when(
      language == "English" ~ "English",
      language == "Unavailable" ~ "Unavailable",
      TRUE ~ "Non-English"
    ), levels = c("English", "Non-English", "Unavailable"))
  )

vars_needed <- c("src_subject_id", "visit", "visit_id", 
                 "severity", "severity_pred", "frequency", "frequency_pred", 
                 "sex", "age_band", "race", "site", "language_group")

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

# Bootstrap ALL Continuous Metrics for one grouping variable
boot_continuous_onevar <- function(data, group_var, ref_group,
                                   n_boot   = N_ITERATIONS,
                                   seed     = SEED,
                                   min_vis  = 10) { 
  
  gsym <- sym(group_var)
  
  d0 <- data %>% filter(!is.na(!!gsym)) %>% mutate(.g = as.character(!!gsym))
  available_groups <- unique(d0$.g)
  
  if (!ref_group %in% available_groups) {
    stop(sprintf("Reference group '%s' not found.", ref_group))
  }
  
  denom <- d0 %>% group_by(.g) %>% summarise(n_visits = n_distinct(visit_id), .groups = "drop")
  keep_groups <- denom %>% filter(n_visits >= min_vis) %>% pull(.g)
  
  if (!ref_group %in% keep_groups) {
    stop(sprintf("Reference group '%s' has fewer than %d visits.", ref_group, min_vis))
  }
  
  ref_n <- denom %>% filter(.g == ref_group) %>% pull(n_visits)
  d0 <- d0 %>% filter(.g %in% keep_groups)
  comp_groups <- sort(setdiff(unique(d0$.g), ref_group))
  
  if (length(comp_groups) == 0) return(tibble())
  
  all_vis_ids <- unique(d0$visit_id)
  
  cat(sprintf("\n[%s] ref='%s' | compare: %s | visits: %d | boots: %d\n",
              group_var, ref_group, paste(comp_groups, collapse = ", "), length(all_vis_ids), n_boot))
  
  # --- RUN WINDOWS-SAFE PARALLEL LOOP WITH PROGRESS BAR ---
  # 'with_progress' tells R to listen for progress updates from the background workers
  boot_results <- with_progress({
    p <- progressor(steps = n_boot)
    
    # future_lapply handles distributing the tasks to the background CPU cores
    future_lapply(seq_len(n_boot), function(b) {
      # Load required libraries inside the worker thread
      library(dplyr)
      library(tibble)
      library(irr)
      
      boot_df <- cluster_boot_sample(d0, all_vis_ids)
      
      metrics_by_g <- boot_df %>%
        group_by(.g) %>%
        group_modify(~ {
          m_sev <- calc_continuous_metrics(.x$severity, .x$severity_pred)
          m_fre <- calc_continuous_metrics(.x$frequency, .x$frequency_pred)
          
          tibble(
            sev_icc = m_sev$icc, sev_mae = m_sev$mae, sev_r = m_sev$pearson,
            fre_icc = m_fre$icc, fre_mae = m_fre$mae, fre_r = m_fre$pearson
          )
        }) %>%
        ungroup()
      
      ref_row <- metrics_by_g %>% filter(.g == ref_group)
      
      res <- list()
      if (nrow(ref_row) == 1) {
        for (g in comp_groups) {
          g_row <- metrics_by_g %>% filter(.g == g)
          if (nrow(g_row) == 1) {
            res[[g]] <- c(
              sev_icc = g_row$sev_icc - ref_row$sev_icc,
              sev_mae = g_row$sev_mae - ref_row$sev_mae,
              sev_r   = g_row$sev_r   - ref_row$sev_r,
              fre_icc = g_row$fre_icc - ref_row$fre_icc,
              fre_mae = g_row$fre_mae - ref_row$fre_mae,
              fre_r   = g_row$fre_r   - ref_row$fre_r
            )
          }
        }
      }
      
      p() # Tick the progress bar forward
      return(res)
      
    }, future.seed = seed) # future.seed guarantees reproducible random sampling across parallel workers
  })
  
  # Extract bootstrap results into matrices for summarization
  boot_sev_icc <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  boot_sev_mae <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  boot_sev_r   <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  
  boot_fre_icc <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  boot_fre_mae <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  boot_fre_r   <- matrix(NA_real_, nrow = n_boot, ncol = length(comp_groups), dimnames = list(NULL, comp_groups))
  
  for (b in seq_len(n_boot)) {
    res_b <- boot_results[[b]]
    if (!is.null(res_b)) {
      for (g in comp_groups) {
        if (!is.null(res_b[[g]])) {
          boot_sev_icc[b, g] <- res_b[[g]]["sev_icc"]
          boot_sev_mae[b, g] <- res_b[[g]]["sev_mae"]
          boot_sev_r[b, g]   <- res_b[[g]]["sev_r"]
          boot_fre_icc[b, g] <- res_b[[g]]["fre_icc"]
          boot_fre_mae[b, g] <- res_b[[g]]["fre_mae"]
          boot_fre_r[b, g]   <- res_b[[g]]["fre_r"]
        }
      }
    }
  }
  
  # Summarize bootstrap results into a tidy data frame with estimates, CIs, and empirical p-values
  summarise_metric <- function(mat, metric_name) {
    lapply(comp_groups, function(g) {
      vals <- mat[, g]
      vals <- vals[!is.na(vals)]
      
      if(length(vals) == 0) return(tibble())
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
    summarise_metric(boot_sev_icc, "Severity: \u0394 ICC"),
    summarise_metric(boot_sev_mae, "Severity: \u0394 MAE"),
    summarise_metric(boot_sev_r,   "Severity: \u0394 Pearson's r"),
    summarise_metric(boot_fre_icc, "Frequency: \u0394 ICC"),
    summarise_metric(boot_fre_mae, "Frequency: \u0394 MAE"),
    summarise_metric(boot_fre_r,   "Frequency: \u0394 Pearson's r")
  ) %>%
    left_join(denom, by = ".g") %>%
    mutate(group_var = group_var, ref = ref_group, ref_n_visits = ref_n)
}


# Setup Parallel Plan and Run Bootstrap
# multisession natively spins up background R sessions on Windows
n_cores <- max(1, future::availableCores() - 1)
plan(multisession, workers = n_cores)
cat(sprintf("\n--- Starting Parallel Processing on %d Cores ---\n", n_cores))

d_sex  <- boot_continuous_onevar(df, "sex",            ref_group = "F")
d_age  <- boot_continuous_onevar(df, "age_band",       ref_group = "18-25")
d_race <- boot_continuous_onevar(df, "race",           ref_group = "White")
d_site <- boot_continuous_onevar(df, "site",           ref_group = "Melbourne") 
d_lang <- boot_continuous_onevar(df, "language_group", ref_group = "English") 

# Shut down background workers when done
plan(sequential) 

# Combine and prepare plotting table
res_df <- bind_rows(
  d_sex  %>% mutate(panel = "Sex"),
  d_age  %>% mutate(panel = "Age band"),
  d_race %>% mutate(panel = "Race"),
  d_site %>% mutate(panel = "Site"),
  d_lang %>% mutate(panel = "Language")
) %>%
  filter(!is.na(estimate)) %>%
  mutate(
    fdr = p.adjust(p_val, method = "BH"),
    fdr_sig = factor(if_else(fdr < 0.05, "Significant (FDR < 0.05)", "Not Significant"), 
                     levels = c("Significant (FDR < 0.05)", "Not Significant")),
    
    panel = factor(panel, levels = c("Sex", "Age band", "Race", "Site", "Language")),
    metric = factor(metric, levels = c("Severity: \u0394 ICC", "Severity: \u0394 MAE", "Severity: \u0394 Pearson's r",
                                       "Frequency: \u0394 ICC", "Frequency: \u0394 MAE", "Frequency: \u0394 Pearson's r")),
    label = paste0(.g, " vs ", ref, "\n(n=", n_visits, ", ref_n=", ref_n_visits, ")")
  ) %>%
  group_by(panel, metric) %>%
  mutate(label = fct_reorder(label, estimate)) %>%
  ungroup()

# Export ALL metrics to CSV (Including Pearson's r)
csv_out <- res_df %>%
  select(panel, group_var, 
         ref_group = ref, reference_n_visits = ref_n_visits, 
         comparison_group = .g, comparison_n_visits = n_visits, 
         metric, estimate, lower_ci_95 = lo, upper_ci_95 = hi, 
         empirical_p_val = p_val, fdr) %>%
  arrange(panel, metric, estimate)

write_csv(csv_out, "bootstrap_continuous_metrics_subgroups_parallel.csv")
cat("\nSaved ALL metrics to: bootstrap_continuous_metrics_subgroups_parallel.csv\n")

# Multi-Panel Forest Plot (EXCLUDING Pearson's r)
plot_df <- res_df %>% 
  filter(!grepl("Pearson's r", metric)) %>%
  mutate(metric = fct_drop(metric))

p_forest <- ggplot(plot_df, aes(x = estimate, y = label, colour = panel, shape = fdr_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_brewer(palette = "Dark2") +
  scale_shape_manual(values = c("Significant (FDR < 0.05)" = 16, "Not Significant" = 1), name = "Significance") +
  facet_grid(panel ~ metric, scales = "free", space = "free_y") +
  labs(
    title    = "Visit-level continuous fairness across demographic subgroups",
    x        = "Difference relative to reference group (95% cluster-bootstrap CI)",
    y        = NULL,
    caption  = sprintf(
      "Cluster bootstrap by visit (n = %d iterations).\nNote: Lower MAE implies lower error. Higher ICC implies higher agreement.\nSignificance adjusted using Benjamini-Hochberg FDR.",
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
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold")
  ) +
  guides(colour = "none") 

print(p_forest)

ggsave("forest_continuous_metrics_grid_parallel.png", p_forest, width = 14, height = 11, dpi = 300)
cat("\nSaved plot to: forest_continuous_metrics_grid_parallel.png\n")