# LLM Validation Analysis Scripts

This directory contains R scripts for comprehensive validation of Large Language Model (LLM) performance on clinical psychiatric assessments, specifically comparing LLM predictions to clinician scores on the Comprehensive Assessment of At-Risk Mental States (CAARMS).

The scripts are mainly designed to be run standalone and reference the data file contained in the data subdirectory of the project.

For the scripts that run with boostrapping, it is advised to run with a low number of boostrap simulations to start off with (say 10) to ensure that the script works.  This can then be upped to the required number of bootstraps to ensure a stable result.

## Overview

The project validates LLM performance across multiple dimensions:
- **Agreement & Reliability**: How well LLM predictions match clinician scores
- **Accuracy**: Balanced accuracy of LLM classifications
- **Fairness & Equity**: Whether performance disparities exist across demographic groups
- **Robustness**: Sensitivity analyses and cross-validation approaches

---

## Script Descriptions

### Core Validation Analyses

#### **LLM Validation bland altman.R**
Performs Bland-Altman agreement analysis comparing LLM vs. clinician continuous scores (severity and frequency ratings). Produces Bland-Altman plots showing bias and limits of agreement, useful for assessing systematic differences and individual-level variability between raters.

#### **LLM Validation continuous scoring metrics.R**
Calculates continuous scoring agreement metrics (intraclass correlation, Pearson correlation, RMSE) between LLM and clinician scores using bootstrap resampling at the visit level. Generates 95% confidence intervals for stable estimates across 10,000 bootstrap iterations.

#### **LLM Validation classification errors.R**
Performs domain-level error analysis identifying which CAARMS domains contribute most to visit-level misclassifications. Maps classification errors to domain-specific drivers to understand root causes of LLM prediction failures.

#### **LLM Validation domain agreement.R**
Assesses agreement between LLM and clinician domain-level classifications. Provides domain-specific balanced accuracy metrics and domain-by-domain reliability statistics to identify which symptom domains the LLM handles well or poorly.

### Fairness & Equity Assessment

#### **LLM Validation propensity score matching.R**
Implements nearest-neighbor propensity score matching (PSM) fairness assessment. Matches participants on baseline characteristics before comparing LLM performance across demographic groups, reducing confounding and estimating causal fairness disparities.

#### **LLM Validation domain fairness significant findings.R**
Conducts domain-level fairness assessment with FDR correction for multiple testing. Outputs only statistically significant findings (p_fdr < 0.05) examining whether LLM performance varies across sex, age, language, site, and race.

#### **LLM Validation evalues psm.R**
Calculates E-values from propensity score matching results. Provides sensitivity analysis to quantify how much confounding would be required to explain away observed fairness disparities, assessing robustness of fairness findings.

### Subgroup Analyses

#### **LLM validation subgroup deltas.R**
Compares LLM performance (balanced accuracy) between subgroups (e.g., females vs. males, English vs. non-English speakers). Estimates delta (difference) statistics with bootstrap confidence intervals to quantify performance disparities. Uses parallelization for efficiency.

#### **LLM validation subgroup deltas - continuous - parallel.R**
Parallelized version calculating continuous scoring metric deltas (e.g., correlation differences) between demographic subgroups. Uses the `future` ecosystem for robust cross-platform parallel processing.

#### **LLM Validation loo multirace fpr.R**
Leave-One-Out (LOO) cross-validation analyzing false positive rates (FPR) across racial groups. Assesses how model generalizes when trained on all races except one, testing robustness across diverse populations.

#### **LLM Validation loo multirace fpr - parallel.R**
Parallelized version of LOO cross-validation with multi-race FPR analysis using base R parallelization and bulletproof progress bars for large-scale computations.

### Aggregated Reporting & Visualizations

#### **LLM Validation with All_Questions data.R**
Generates domain-level balanced accuracy heatmaps stratified by demographics (sex, age band, language, site, race). Each heatmap cell shows mean balanced accuracy with 95% CI and sample size for subgroup-domain combinations. Uses CAARMS CHR criterion (APS2a/APS2b/BLIPS rules).

#### **LLM Validation overall symptom breakdown - BA.R**
Produces Bland-Altman agreement visualizations broken down by overall symptom severity categories. Examines agreement patterns across symptom severity strata and generates publication-quality plots.

#### **LLM Validation Generic with Plot Types and CI.R**
Generic validation script with flexible plot types and confidence interval generation. Provides configurable subgroup analysis by sex, age band, and language with customizable predictor categories for broader applicability.

#### **LLM Validation WPTU analysis.R**
Site-specific analysis for Western Psychiatric Treatment Unit (WPTU). Analyzes LLM performance trends using multi-race FPR calculations and compares ME (likely Maine Medical Center) vs. non-ME sites.

### Data Quality & Integrity

#### **LLM Validation validation check.R**
Core validation checks ensuring data integrity before analysis:
- Loads and preprocesses phenotype data
- Validates visit-level aggregation (expected 547 visits from 331 unique subjects)
- Applies standard factor encoding and derived variables
- Serves as prerequisite for other analyses

#### **LLM Validation subgroup  sample size audit.R**
Audits sample sizes for subgroup analyses. Identifies which demographic-domain combinations meet minimum subject thresholds and ensures adequate power for statistical inference.

#### **LLM Validation WPTU analysis.R** (also includes data inspection)
Includes Shapiro-Wilk normality testing and baseline severity distribution comparisons between included vs. excluded participants, identifying potential selection bias.

---

## Data Files

All scripts expect data at:
- **Primary dataset**: `data/complete_consolidated_sheets.xlsx` 
  - Sheet: `"No Missing and No Partial"` (main analytic sample)
  - Also uses `"No Partial"` sheet for baseline comparisons

### Required Data Columns
- `src_subject_id`: Subject identifier
- `visit`: Visit number/identifier
- `question` (Q1-Q16): CAARMS question IDs
- `severity`, `frequency`: Clinician ratings (numeric)
- `severity_pred`, `frequency_pred`: LLM predictions (numeric)
- `sex`, `age_band`, `language`, `race`: Demographic variables
- `interview_date`: Visit timestamp

---

## Helper Functions & Dependencies

### Helper Scripts (sourced by analyses)
- `bootstrap_helpers.R`: Bootstrap resampling utilities and statistical functions
- `visit_level_helpers.R`: Aggregation from domain-level to visit-level classifications

### Key R Libraries
- **Data manipulation**: `dplyr`, `tidyr`, `readxl`
- **Statistical inference**: `irr`, `boot`, `MatchIt`, `Metrics`
- **Visualization**: `ggplot2`, `forcats`
- **Parallelization**: `future`, `furrr`, `future.apply`, `pbapply`, `parallel`
- **Reporting**: `knitr`, `kableExtra`, `readr`

---

## CAARMS Criteria Reference

Scripts use standardized CAARMS thresholds for At-Risk Mental State (APS) classification:

| Criterion | Severity Range | Frequency Range |
|-----------|-----------------|-----------------|
| **APS2a** | 3-5 | 3-6 |
| **APS2b** | 6 | 3 |
| **BLIPS** | 6 | 4-6 |

CHR (Clinical High Risk) designation: Meets at least one APS criterion (APS2a, APS2b, or BLIPS) or BLIPS criterion across analyzed domains.

---

## Typical Analysis Workflow

1. **Start**: `LLM Validation validation check.R` - Ensure data integrity
2. **Agreement**: `LLM Validation continuous scoring metrics.R` + `LLM Validation bland altman.R`
3. **Classification**: `LLM Validation domain agreement.R` + `LLM Validation classification errors.R`
4. **Fairness**: `LLM Validation propensity score matching.R` → `LLM Validation evalues psm.R`
5. **Subgroup Disparities**: `LLM validation subgroup deltas.R` + `LLM Validation loo multirace fpr.R`
6. **Reporting**: `LLM Validation with All_Questions data.R` + `LLM Validation overall symptom breakdown - BA.R`
7. **Audit**: `LLM Validation subgroup  sample size audit.R`

---

## Output Notes

Most scripts generate:
- **Tables/CSVs**: Summary statistics and fairness metrics
- **Plots**: Publication-quality visualizations (Bland-Altman, heatmaps, forest plots)
- **Console output**: Progress tracking and validation messages

Bootstrap-based analyses use consistent `SEED = 13` and typically `N_BOOT = 10000` iterations for stable confidence intervals.

---

