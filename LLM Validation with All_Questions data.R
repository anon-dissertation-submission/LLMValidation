# =============================================================================
# Domain-level Balanced Accuracy heatmaps — visit-level bootstrapped
#   Heatmap 1: Sex (F vs M)    × 15 symptom domains
#   Heatmap 2: Age band        × 15 symptom domains
#   Heatmap 3: Language        × 15 symptom domains
#   Heatmap 4: Site            × 15 symptom domains
#   Heatmap 5: Race            × 15 symptom domains
#
# Fill = bootstrap MEAN balanced accuracy (resampled by visit)
# CHR criterion per domain: full CAARMS rule (APS2a / APS2b / BLIPS)
# Labels: Mean BA [95% CI] (n_subjects)
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(readr)
library(tibble)

# Settings
DATA_FILE    <- "data/complete_consolidated_sheets.xlsx"
DATA_SHEET   <- "No Missing and No Partial"

N_BOOT       <- 10000   # Set to 1000 for final production
SEED         <- 13
MIN_SUBJ     <- 11     # Minimum unique subjects per cell (drops groups where n <= 10)

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

# Domain map
domain_map <- tibble::tribble(
  ~q_num, ~domain,
  1L,  "Unusual thoughts & experiences",
  2L,  "Suspiciousness",
  3L,  "Unusual somatic ideas",
  4L,  "Ideas of guilt",
  5L,  "Jealous ideas",
  6L,  "Unusual religious ideas",
  7L,  "Erotomanic ideas",
  8L,  "Grandiosity",
  9L,  "Auditory perceptual abnormalities",
  10L, "Visual perceptual abnormalities",
  11L, "Olfactory perceptual abnormalities",
  12L, "Gustatory perceptual abnormalities",
  13L, "Tactile perceptual abnormalities",
  14L, "Somatic perceptual abnormalities",
  15L, "Disorganised communication / expression"
)

# Helpers
is_chr_domain <- function(sev, freq) {
  (sev %in% CAARMS_SEVERITY_LOW_APS2a:CAARMS_SEVERITY_HIGH_APS2a &
     freq %in% CAARMS_FREQUENCY_LOW_APS2a:CAARMS_FREQUENCY_HIGH_APS2a) |
    (sev == CAARMS_SEVERITY_APS2b  & freq == CAARMS_FREQUENCY_APS2b) |
    (sev == CAARMS_SEVERITY_BLIPS  & freq %in% CAARMS_FREQUENCY_LOW_BLIPS:CAARMS_FREQUENCY_HIGH_BLIPS)
}

calc_ba <- function(truth, pred) {
  truth <- as.integer(truth); pred <- as.integer(pred)
  n_pos <- sum(truth == 1, na.rm = TRUE)
  n_neg <- sum(truth == 0, na.rm = TRUE)
  tpr   <- if (n_pos > 0) sum(pred == 1 & truth == 1, na.rm = TRUE) / n_pos else NA_real_
  tnr   <- if (n_neg > 0) sum(pred == 0 & truth == 0, na.rm = TRUE) / n_neg else NA_real_
  mean(c(tpr, tnr), na.rm = TRUE)
}

visit_boot_sample <- function(data, visit_ids) {
  sampled <- sample(visit_ids, size = length(visit_ids), replace = TRUE)
  tbl     <- tibble(visit_id = sampled, .rep = seq_along(sampled))
  data %>%
    inner_join(tbl, by = "visit_id", relationship = "many-to-many") %>%
    select(-.rep)
}

# Load + prep data
raw <- read_excel(DATA_FILE, sheet = DATA_SHEET) %>%
  mutate(
    src_subject_id = as.character(src_subject_id),
    visit = as.character(visit),
    visit_id = paste0(src_subject_id, "_", visit),
    sex      = factor(sex),
    age_band = factor(age_band),
    lang_group = case_when(
      is.na(language) ~ "Unavailable",
      tolower(trimws(language)) == "english" ~ "English",
      tolower(trimws(language)) == "unavailable" ~ "Unavailable",
      TRUE ~ "Non-English"
    ),
    lang_group = factor(lang_group, levels = c("English", "Non-English", "Unavailable")),
    site_group = case_when(
      is.na(site) ~ NA_character_,
      toupper(trimws(site)) == "ME" ~ "ME",
      TRUE ~ "Not-ME"
    ),
    site_group = factor(site_group, levels = c("ME", "Not-ME")),
    race_group = factor(trimws(race), levels = c("White", "Asian", "Black or African American", 
                                                 "More than one race", "American Indian/Alaska Native", 
                                                 "Unknown or not reported"))
  )

domain_rows <- raw %>%
  mutate(q_num = readr::parse_number(as.character(question))) %>%
  left_join(domain_map, by = "q_num") %>%
  filter(
    !is.na(domain),
    !is.na(severity),      !is.na(frequency),
    !is.na(severity_pred), !is.na(frequency_pred)
  ) %>%
  mutate(
    chr_true = as.integer(is_chr_domain(severity,      frequency)),
    chr_pred = as.integer(is_chr_domain(severity_pred, frequency_pred))
  )

# Core bootstrap function
boot_ba_by_domain <- function(data, group_var,
                              n_boot   = N_BOOT,
                              seed     = SEED,
                              min_subj = MIN_SUBJ) {
  set.seed(seed)
  gsym <- rlang::sym(group_var)
  
  d0 <- data %>%
    filter(!is.na(!!gsym)) %>%
    mutate(.grp = as.character(!!gsym))
  
  groups  <- sort(unique(d0$.grp))
  domains <- sort(unique(d0$q_num))
  
  subj_per_cell <- d0 %>%
    group_by(.grp, q_num) %>%
    summarise(n_subj = n_distinct(src_subject_id), .groups = "drop")
  
  all_visit_ids <- unique(d0$visit_id)
  
  boot_store <- lapply(groups, function(g) {
    matrix(NA_real_, nrow = n_boot, ncol = length(domains),
           dimnames = list(NULL, as.character(domains)))
  })
  names(boot_store) <- groups
  
  for (b in seq_len(n_boot)) {
    boot_df <- visit_boot_sample(d0, all_visit_ids)
    for (g in groups) {
      g_df <- boot_df %>% filter(.grp == g)
      if (nrow(g_df) == 0) next
      for (qn in domains) {
        q_df <- g_df %>% filter(q_num == qn)
        if (nrow(q_df) < 2) next
        ba <- calc_ba(q_df$chr_true, q_df$chr_pred)
        boot_store[[g]][b, as.character(qn)] <- ba
      }
    }
  }
  
  lapply(groups, function(g) {
    mat <- boot_store[[g]]
    lapply(as.character(domains), function(qn) {
      vals <- mat[, qn]; vals_ok <- vals[!is.na(vals)]
      n_s  <- subj_per_cell %>% filter(.grp == g, q_num == as.integer(qn)) %>% pull(n_subj)
      tibble(group_var = group_var, group = g, q_num = as.integer(qn), n_subj = if(length(n_s)>0) n_s else 0,
             n_valid = length(vals_ok), ba_mean = mean(vals_ok), ba_lo = quantile(vals_ok, 0.025), ba_hi = quantile(vals_ok, 0.975))
    }) %>% bind_rows()
  }) %>% bind_rows() %>% left_join(domain_map, by = "q_num")
}

# Run bootstrap
cat("Running Sex...\n")
ba_sex  <- boot_ba_by_domain(domain_rows, "sex")
cat("Running Age...\n")
ba_age  <- boot_ba_by_domain(domain_rows %>% filter(!is.na(age_band)), "age_band")
cat("Running Language...\n")
ba_lang <- boot_ba_by_domain(domain_rows, "lang_group")
cat("Running Site...\n")
ba_site <- boot_ba_by_domain(domain_rows %>% filter(!is.na(site_group)), "site_group")
cat("Running Race...\n")
ba_race <- boot_ba_by_domain(domain_rows %>% filter(!is.na(race_group)), "race_group")

# Labels + Processing
domain_order <- domain_rows %>%
  group_by(q_num) %>%
  summarise(n_total = n_distinct(src_subject_id), .groups = "drop") %>%
  left_join(domain_map, by = "q_num") %>%
  arrange(n_total) %>%
  mutate(domain_lab = domain)

process_for_plot <- function(df, group_levels = NULL) {
  # 1. Identify groups that have at least MIN_SUBJ (11) subjects in at least one cell
  valid_groups <- df %>%
    group_by(group) %>%
    summarise(max_n = max(n_subj, na.rm = TRUE), .groups = "drop") %>%
    filter(max_n >= MIN_SUBJ) %>%
    pull(group)
  
  df %>%
    # 2. Filter out groups entirely if they fall below the threshold
    filter(group %in% valid_groups) %>%
    mutate(domain_lab = factor(domain, levels = domain_order$domain_lab)) %>%
    {if (!is.null(group_levels)) {
      # Keep only the levels of groups that passed the filter, so empty columns drop off
      mutate(., group = factor(group, levels = intersect(group_levels, valid_groups))) 
    } else {
      mutate(., group = factor(group))
    }} %>%
    mutate(
      ba_plot    = if_else(n_subj < MIN_SUBJ | n_valid < (N_BOOT * 0.5), NA_real_, ba_mean),
      tile_label = if_else(is.na(ba_plot), "—", sprintf("%.2f\n[%.2f, %.2f]\n(n=%d)", ba_plot, ba_lo, ba_hi, n_subj)),
      txt_col    = if_else(is.na(ba_plot), "grey60", if_else(abs(ba_plot - 0.5) > 0.25, "white", "black"))
    )
}

ba_sex  <- process_for_plot(ba_sex)
ba_age  <- process_for_plot(ba_age, c("<18", "18-25", ">25"))
ba_lang <- process_for_plot(ba_lang, c("English", "Non-English", "Unavailable"))
ba_site <- process_for_plot(ba_site, c("ME", "Not-ME"))
ba_race <- process_for_plot(ba_race, c("White", "Asian", "Black or African American", 
                                       "More than one race", "American Indian/Alaska Native", 
                                       "Unknown or not reported"))

# Plotting
ba_fill_scale <- scale_fill_gradient2(
  low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", 
  midpoint = 0.50, limits = c(0, 1), name = "BA Mean"
)

heatmap_theme <- theme_minimal(base_size = 12) + 
  theme(
    panel.grid   = element_blank(), 
    axis.text.x  = element_text(face = "bold", size = 12, colour = "black"), 
    axis.text.y  = element_text(face = "bold", size = 11, colour = "black"), 
    plot.title   = element_text(face = "bold", size = 13),
    plot.caption = element_text(size = 8, hjust = 0, colour = "grey50")
  )

heatmap_theme_race <- heatmap_theme + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# Updated title text block
plot_subtitle <- "Tiles show Mean [95% CI] and Sample Size (n) and excluding groups with <= 10 subjects"
boot_caption  <- "Fill = visit-level bootstrap mean BA."
boot_caption  <- sprintf("Fill = visit-level bootstrap mean BA. Labels: Mean [95%% CI] (n=subjects).\nGroups with n <= %d are excluded.", MIN_SUBJ - 1)

p_sex  <- ggplot(ba_sex, aes(x = group, y = domain_lab, fill = ba_plot)) + geom_tile(colour="white") + geom_text(aes(label=tile_label, colour=txt_col), size=3, fontface="bold") + scale_colour_identity() + ba_fill_scale + labs(title="BA by Symptom Domain and Sex", caption=boot_caption, x=NULL, y=NULL) + heatmap_theme
p_age  <- ggplot(ba_age, aes(x = group, y = domain_lab, fill = ba_plot)) + geom_tile(colour="white") + geom_text(aes(label=tile_label, colour=txt_col), size=3, fontface="bold") + scale_colour_identity() + ba_fill_scale + labs(title="BA by Symptom Domain and Age Band", caption=boot_caption, x=NULL, y=NULL) + heatmap_theme
p_lang <- ggplot(ba_lang, aes(x = group, y = domain_lab, fill = ba_plot)) + geom_tile(colour="white") + geom_text(aes(label=tile_label, colour=txt_col), size=3, fontface="bold") + scale_colour_identity() + ba_fill_scale + labs(title="BA by Symptom Domain and Language", caption=boot_caption, x=NULL, y=NULL) + heatmap_theme
p_site <- ggplot(ba_site, aes(x = group, y = domain_lab, fill = ba_plot)) + geom_tile(colour="white") + geom_text(aes(label=tile_label, colour=txt_col), size=3, fontface="bold") + scale_colour_identity() + ba_fill_scale + labs(title="BA by Symptom Domain and Site", caption=boot_caption, x=NULL, y=NULL) + heatmap_theme
p_race <- ggplot(ba_race, aes(x = group, y = domain_lab, fill = ba_plot)) + geom_tile(colour="white") + geom_text(aes(label=tile_label, colour=txt_col), size=3, fontface="bold") + scale_colour_identity() + ba_fill_scale + labs(title="BA by Symptom Domain and Race", caption=boot_caption, x=NULL, y=NULL) + heatmap_theme_race

# Save

cat("\nExporting BA data to CSVs...\n")

readr::write_csv(ba_sex,  "ba_summary_sex.csv")
readr::write_csv(ba_age,  "ba_summary_age.csv")
readr::write_csv(ba_lang,  "ba_summary_language.csv")
readr::write_csv(ba_site,  "ba_summary_site.csv")
readr::write_csv(ba_race,  "ba_summary_race.csv")

cat("Exports complete.\n\n")

cat("\nGenerating heatmaps...\n")
ggsave("heatmap_BA_sex.png", p_sex, width = 7, height = 9)
ggsave("heatmap_BA_age.png", p_age, width = 8, height = 9)
ggsave("heatmap_BA_language.png", p_lang, width = 8, height = 9)
ggsave("heatmap_BA_site.png", p_site, width = 7, height = 9)
ggsave("heatmap_BA_race.png", p_race, width = 12, height = 9)
cat("Heatmaps complete.\n\n")

