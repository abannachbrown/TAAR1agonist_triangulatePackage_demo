#### TAAR-1 Triangulation Example

### Load packages. Additional plotting helpers are sourced from the local
### Functions/ folder below, with the used functions documented there.
library(triangulate)
library(dplyr)
library(readr)
library(metafor)
library(meta)
library(ggplot2)
suppressPackageStartupMessages(library(stringr))

###############################################################################
### Source local plotting helper functions
###
### Functions/helpers.R
### - clean_data(): converts risk-of-bias judgements, bias types, and absolute
###   bias directions to the one-letter plotting codes used by rob_direction().
###
source("Functions/helpers.R")

### Functions/helpers-metafor.R
### - .pval(): formats meta-analysis p-values in model summary labels.
### - .fcf() and .psort(): format and order estimate / confidence-interval text.
### - annotate_poly(): prints pooled estimate [95% CI] text for subgroup summary
###   polygons drawn by metafor::addpoly().
###
source("Functions/helpers-metafor.R")

### Functions/BiasCorrect_robviz.R - custom additions to "triangulate" package 
### - rob_direction(): local amended plotting function used at the end of this
###   script. It draws the original forest plot, risk-of-bias traffic-light
###   columns, and bias-corrected forest plot together. This local version has
###   been amended from the package internals to show D1-D10, include Animal and
###   RCT subgroups, keep CIs away from the risk-of-bias grid, and use
###   ASCII-safe labels for reliable PDF output.
### - atransf argument: the original package-style plotting code used exp(),
###   which is appropriate for log-ratio outcomes but wrong for SMDs. The local
###   function now defaults to identity and is called explicitly with identity
###   below, so original yi and bias-adjusted yi_adj are displayed on the SMD
###   scale without exponentiation.
source("Functions/BiasCorrect_robviz.R")

### Functions/animal-helpers.R
### - subgroup_analysis(): recreates the original animal drug-level multivariate
###   meta-analysis used for the 17 visible animal rows.
### - forest_subgroup(): optional ggplot forest plot for checking the original
###   animal model outside the triangulation plot.
### - create_formula(): creates the nested random-effects formula used in the
###   original 188-effect animal rma.mv() models.
### from LSR3-animals #### from https://github.com/galenos-project/LSR3_taar1_A/blob/main/LSR3_animal_analysis_u1.html  
source("Functions/animal-helpers.R")



#### load original animal df to recreate the multi-variate model from 188 experiments
#### from https://github.com/galenos-project/LSR3_taar1_A/blob/main/LSR3_animal_analysis_u1.html 
df <- read_csv(
  "animal_df_full.csv"
)


#Figure 2.1.4.7 - Subgroup analysis of TAAR1 agonist vs control on locomotor activity by intervention administered
SMD_S_LMA_Drug <- subgroup_analysis(df, "TvC", "Locomotor activity", "DrugName", 0.5)
# Optional check against the original ggplot forest:
# forest_subgroup(SMD_S_LMA_Drug$plotdata, "DrugName", "Locomotor Activity", "Drug")

##################
### Load animal studies

# taar1_drug_new <- read_csv("taar1_drug_merged.csv")
# head(taar1_drug_new)
     
taar1_animal <- read_csv("taar1_drug_merged.csv") %>%
  mutate(
    # tri_prep_data() aggregates by result_id, so result_id must identify one
    # effect estimate. The new animal CSV has one row per drug-level estimate.
    result_id = DrugName
  ) %>%
  rename(study = DrugName) %>%
  mutate(
    type = "Animal",
    # Use the same judgement vocabulary for animal and human studies.
    # SYRCLE-style coding: Yes = low risk, No = high risk, Unclear = moderate.
    overall = "Moderate"
  )

# The 17-row animal CSV contains RoB judgement columns and BiasType columns.
# If BiasDirection columns are absent, add the planned domain-level directions
# used by triangulate before renaming to d{domain}{j/t/d}.
animal_bias_direction <- tibble::tibble(
  Domain = c(
    "D1 Allocation sequence",
    "D2 Baseline similarity",
    "D3 Allocation concealment",
    "D4 Random housing",
    "D5 Caregivers blinded",
    "D6 Random selection for outcome assessment",
    "D7 Blinded outcome assessor",
    "D8 Incomplete data reporting addressed",
    "D9 Free from selective outcome reporting",
    "D10 Free from other risks of bias"
  ),
  BiasDirection = c(
    "Away from null",
    "Unpredictable",
    "Away from null",
    "Unpredictable",
    "Away from null",
    "Unpredictable",
    "Away from null",
    "Unpredictable",
    "Away from null",
    "Unpredictable"
  )
)

for (i in seq_len(nrow(animal_bias_direction))) {
  direction_col <- paste0(animal_bias_direction$Domain[i], " BiasDirection")
  if (!direction_col %in% names(taar1_animal)) {
    taar1_animal[[direction_col]] <- animal_bias_direction$BiasDirection[i]
  }
}

# Convert verbose animal domain columns from the CSV into triangulate's
# d{domain}{j/t/d} format, e.g. d1j, d1t, d1d.
domain_cols <- grep("^D[0-9]+ ", names(taar1_animal), value = TRUE)

names(taar1_animal)[match(domain_cols, names(taar1_animal))] <- vapply(
  domain_cols,
  function(col) {
    domain_num <- str_match(col, "^D([0-9]+) ")[, 2]
    suffix <- case_when(
      str_detect(col, " BiasType$") ~ "t",
      str_detect(col, " BiasDirection$") ~ "d",
      TRUE ~ "j"
    )
    paste0("d", domain_num, suffix)
  },
  character(1)
)

taar1_animal <- taar1_animal %>%
  mutate(
    across(matches("^d[0-9]+j$"), ~ case_when(
      .x == "Yes" ~ "Low",
      .x == "No" ~ "High",
      .x == "Unclear" ~ "Moderate",
      is.na(.x) ~ "None",
      TRUE ~ .x
    )),
    overall = apply(
      select(., matches("^d[0-9]+j$")),
      1,
      function(x) {
        x <- unique(stats::na.omit(as.character(x)))
        if ("High" %in% x) return("High")
        if ("Moderate" %in% x) return("Moderate")
        if ("Low" %in% x) return("Low")
        "None"
      }
    )
  )


##################
### Load human RCTs

### human_taar1.csv contains the meta-analysis effect estimate columns
### directly:
### - yi = study-level SMD
### - se = study-level standard error
### - vi = study-level variance, equal to se^2
### It also contains the ROB judgements joined by the shortened study label.
taar1_human_4 <- read_csv("human_taar1.csv") %>%
  mutate(
    result_id = as.character(row_number()),
    type = "RCT"
  ) %>%
  rename(
    d1j = "Bias due to randomization",
    d2j = "Bias due to deviations from intended intervention",
    d3j = "Bias due to missing data",
    d4j = "Bias due to outcome measurement",
    d5j = "Bias due to selected reported results",
    d6j = "Bias due to period and carryover effects",
    overall = "Overall bias"
  ) %>%
  mutate(
    # Human ROB2-style judgements already use Low / Some concerns / High.
    # Convert "Some concerns" to "Moderate" so it matches the plotting and
    # custom-prior scale used for the animal studies.
    across(c(d1j, d2j, d3j, d4j, d5j, d6j, overall), ~ case_when(
      .x == "Some concerns" ~ "Moderate",
      is.na(.x) ~ "None",
      TRUE ~ as.character(.x)
    )),

    # Bias types. "prop" means a proportional bias prior; "None" means no
    # adjustment for that absent/not-applicable domain.
    d1t = "prop",
    d2t = "prop",
    d3t = "prop",
    d4t = "prop",
    d5t = "prop",
    d6t = "None",

    # Relative directions of bias. tri_absolute_direction() converts these to
    # absolute left/right directions using yi and the bias type.
    d1d = "Unpredictable",
    d2d = "Towards null",
    d3d = "Away from null",
    d4d = "Away from null",
    d5d = "Away from null",
    d6d = "None",

    # Add empty D7-D10 columns so the amended D1-D10 plotting function can draw
    # human and animal studies in the same figure.
    d7j = "None",  d7t = "None",  d7d = "None",
    d8j = "None",  d8t = "None",  d8d = "None",
    d9j = "None",  d9t = "None",  d9d = "None",
    d10j = "None", d10t = "None", d10d = "None"
  )

if (any(is.na(taar1_human_4$yi)) || any(is.na(taar1_human_4$vi))) {
  stop("human_taar1.csv must contain non-missing yi and vi columns.")
}

if ("se" %in% names(taar1_human_4) &&
    any(abs(taar1_human_4$vi - taar1_human_4$se^2) > 1e-6, na.rm = TRUE)) {
  stop("human_taar1.csv vi must equal se^2.")
}


##################
### Check the long-format data expected by triangulate

taar1_animal_long <- taar1_animal %>%
  tri_to_long()

taar1_human_long <- taar1_human_4 %>%
  tri_to_long()

triangulate::tri_dat_check(taar1_animal_long, mode = "full")
triangulate::tri_dat_check(taar1_human_long, mode = "full")

# These objects are retained so you can inspect how relative directions have
# been converted before priors are joined.
taar1_animal_abs <- taar1_animal_long %>%
  tri_absolute_direction()

taar1_human_abs <- taar1_human_long %>%
  tri_absolute_direction()


##################
### Indirectness data

### There are no bespoke indirectness assessments in the current TAAR1 files.
### We therefore create explicit "None" indirectness rows for D1-D3 for each
### effect estimate, then still run the standard triangulate functions so the
### adjustment pipeline is transparent.

taar1_animal_ind <- taar1_animal %>%
  select(result_id, study, type, yi, vi) %>%
  mutate(
    d1j = "None", d1t = "None", d1d = "None",
    d2j = "None", d2t = "None", d2d = "None",
    d3j = "None", d3t = "None", d3d = "None"
  ) %>%
  tri_to_long() %>%
  tri_absolute_direction() %>%
  tri_append_indirect(triangulate::dat_ind_values)

taar1_human_ind <- taar1_human_4 %>%
  select(result_id, study, type, yi, vi) %>%
  mutate(
    d1j = "None", d1t = "None", d1d = "None",
    d2j = "None", d2t = "None", d2d = "None",
    d3j = "None", d3t = "None", d3d = "None"
  ) %>%
  tri_to_long() %>%
  tri_absolute_direction() %>%
  tri_append_indirect(triangulate::dat_ind_values)


##################
### Custom bias priors for animal studies

### Animal studies use the harmonised Low / Moderate / High scale created above.
### High risk is anchored to the animal pooled SMD half-width:
### (1.38 - 0.84) / 2 = 0.27.
### The proportional high-risk prior is 0.27 / 1.11 = about 0.24.
###
### Additive high risk: N(0.27, 0.018)
### SD ~= 0.134, so about 95% of mass lies roughly from 0.01 to 0.53.
###
### Proportional high risk: N(0.24, 0.014)
### SD ~= 0.118, so about 95% of mass lies roughly from 0.01 to 0.47.
###
### Practical scaling rule:
### low      = 0.25 * high
### moderate = 0.50 * high
### high     = 1.00 * high
### critical = 2.00 * high
###
### The same multiplier is applied to both means and variances, consistent with
### the triangulate convention used for more severe bias categories.

custom_bias_priors_animal <- tibble::tibble(
  domain = "all",
  j = c("low", "moderate", "high", "critical"),
  # ADDITIVE PRIORS
  # low:      mean = 0.0675  (25% of high risk)
  # moderate: mean = 0.135   (50% of high risk)
  # high:     mean = 0.27    (animal pooled SMD CI half-width)
  # critical: mean = 0.54    (2 * high risk)
  bias_m_add = c(0.0675, 0.135, 0.27, 0.54),
  bias_v_add = c(0.0045, 0.009, 0.018, 0.036),
  # PROPORTIONAL PRIORS
  # low:      mean = 0.06    (25% of high risk)
  # moderate: mean = 0.12    (50% of high risk)
  # high:     mean = 0.24    (0.27 / pooled SMD 1.11)
  # critical: mean = 0.48    (2 * high risk)
  bias_m_prop = c(0.06, 0.12, 0.24, 0.48),
  bias_v_prop = c(0.0035, 0.007, 0.014, 0.028)
)

taar1_animal_prepped <- taar1_animal %>%
  tri_to_long() %>%
  tri_absolute_direction() %>%
  tri_append_bias(custom_bias_priors_animal)


##################
### Custom bias priors for human RCTs

### Human studies use a separate prior scale because the empirical anchor is
### different. High risk is anchored to the human pooled SMD CI half-width:
### (0.34 - (-0.05)) / 2 = 0.195.
### The proportional high-risk prior is capped at 0.50 because the pooled
### effect is close to the null and the raw ratio would be unstable.

custom_bias_priors_human <- tibble::tibble(
  domain = "all",
  j = c("low", "moderate", "high", "critical"),
  # LOW RISK OF BIAS
  # Additive mean = 0.04875 (25% of high risk: 0.195)
  # Represents small residual bias; anchored to CI half-width (0.195) of pooled SMD.
  # Assumes low-risk studies may still have minimal bias, but unlikely to materially shift effect.
  #
  # MODERATE RISK OF BIAS
  # Additive mean = 0.0975 (50% of high risk)
  # Represents moderate bias, roughly half the CI half-width.
  # Could meaningfully shift effect estimates but not dominate them.
  #
  # HIGH RISK OF BIAS
  # Additive mean = 0.195 (CI half-width of pooled SMD: (0.34 - (-0.05)) / 2)
  # Anchors bias magnitude to empirical uncertainty in the meta-analysis.
  # Represents bias large enough to shift estimate across null.
  #
  # CRITICAL / VERY HIGH RISK OF BIAS
  # Additive mean = 0.39 (2 * 0.195)
  # Reflects extreme bias scenarios and allows bias to dominate observed effect.
  bias_m_add = c(0.04875, 0.0975, 0.195, 0.39),
  # Additive variances:
  # low      = 0.005, scaled proportionally from high variance 0.02
  # moderate = 0.01, scaled proportionally
  # high     = 0.02, chosen to reflect moderate uncertainty (SD ~= 0.14)
  # critical = 0.04, doubled from high risk
  # The high-risk variance is consistent with observed heterogeneity (tau ~= 0.157).
  bias_v_add = c(0.005, 0.01, 0.02, 0.04),
  # Proportional means:
  # low      = 0.125 (25% of high risk: 0.50)
  # moderate = 0.25  (50% of high risk)
  # high     = 0.50  (capped rather than raw ratio ~= 1.30)
  # critical = 1.00  (2 * high risk)
  # The cap avoids instability when the pooled effect is near null.
  bias_m_prop = c(0.125, 0.25, 0.50, 1.00),
  # Proportional variances:
  # low      = 0.02, scaled proportionally from high variance 0.08
  # moderate = 0.04, scaled proportionally
  # high     = 0.08, allows substantial uncertainty (SD ~= 0.28)
  # critical = 0.16, doubled from high risk
  bias_v_prop = c(0.02, 0.04, 0.08, 0.16)
)

taar1_human_prepped <- taar1_human_4 %>%
  tri_to_long() %>%
  tri_absolute_direction() %>%
  tri_append_bias(custom_bias_priors_human)


##################
### Calculate adjusted estimates

taar1_animal_final <- tri_prep_data(
  taar1_animal_prepped,
  taar1_animal_ind
)

taar1_human_final <- tri_prep_data(
  taar1_human_prepped,
  taar1_human_ind
)


##################
### Animal multivariate summary models for the triangulation plot

### The triangulation plot shows 17 animal drug-level rows, but the animal
### subgroup summary should use the original 188-effect multivariate model:
### rma.mv() with vcalc() covariance and nested random effects. The block below
### recreates that model for the original panel and then creates an analogous
### bias-adjusted 188-effect model for the adjusted panel.

animal_mv_base <- df %>%
  filter(SortLabel == "TvC") %>%
  filter(outcome_type == "Locomotor activity") %>%
  filter(!is.na(SMDv)) %>%
  filter(!is.na(DrugName)) %>%
  filter(DrugName %in% taar1_animal$study) %>%
  mutate(
    effect_id = row_number(),
    result_id = paste0("animal_mv_", effect_id),
    study = DrugName,
    type = "Animal",
    yi = SMD,
    vi = SMDv
  )

animal_mv_random_formula <- create_formula(
  c("Strain", "StudyId", "ExperimentID_I"),
  animal_mv_base
)

if (is.null(animal_mv_random_formula)) {
  stop("Unable to create animal multivariate random-effects formula.")
}

animal_mv_v <- metafor::vcalc(
  vi = vi,
  cluster = StudyId,
  subgroup = ExperimentID_I,
  obs = effect_id,
  data = animal_mv_base,
  rho = 0.5
)

animal_mv_original_model <- metafor::rma.mv(
  yi = yi,
  V = animal_mv_v,
  random = animal_mv_random_formula,
  data = animal_mv_base,
  method = "REML",
  test = "t",
  dfs = "contain",
  control = list(optimizer = "nlminb")
)

animal_mv_rob <- taar1_animal %>%
  select(study, matches("^d[0-9]+(j|t|d)$"), overall)

animal_mv_for_tri <- animal_mv_base %>%
  select(result_id, study, type, yi, vi, effect_id, StudyId, ExperimentID_I, Strain) %>%
  left_join(animal_mv_rob, by = "study")

if (any(is.na(animal_mv_for_tri$d1j))) {
  stop("Some 188-effect animal rows failed to join to the 17-row drug-level RoB table.")
}

animal_mv_bias <- animal_mv_for_tri %>%
  select(result_id, study, type, yi, vi, matches("^d[0-9]+(j|t|d)$")) %>%
  tri_to_long() %>%
  tri_absolute_direction() %>%
  tri_append_bias(custom_bias_priors_animal)

animal_mv_ind <- animal_mv_for_tri %>%
  select(result_id, study, type, yi, vi) %>%
  mutate(
    d1j = "None", d1t = "None", d1d = "None",
    d2j = "None", d2t = "None", d2d = "None",
    d3j = "None", d3t = "None", d3d = "None"
  ) %>%
  tri_to_long() %>%
  tri_absolute_direction() %>%
  tri_append_indirect(triangulate::dat_ind_values)

animal_mv_final <- tri_prep_data(animal_mv_bias, animal_mv_ind)

animal_mv_adjusted_data <- animal_mv_base %>%
  left_join(
    animal_mv_final %>% select(result_id, yi_adj, vi_adj),
    by = "result_id"
  )

animal_mv_v_adj <- metafor::vcalc(
  vi = vi_adj,
  cluster = StudyId,
  subgroup = ExperimentID_I,
  obs = effect_id,
  data = animal_mv_adjusted_data,
  rho = 0.5
)

animal_mv_adjusted_model <- metafor::rma.mv(
  yi = yi_adj,
  V = animal_mv_v_adj,
  random = animal_mv_random_formula,
  data = animal_mv_adjusted_data,
  method = "REML",
  test = "t",
  dfs = "contain",
  control = list(optimizer = "nlminb")
)

animal_external_original_summary <- tibble::tibble(
  type = "Animal",
  yi = as.numeric(animal_mv_original_model$beta),
  ci.lb = as.numeric(animal_mv_original_model$ci.lb),
  ci.ub = as.numeric(animal_mv_original_model$ci.ub),
  label = "Animal rma.mv model (188 effects)"
)

animal_external_adjusted_summary <- tibble::tibble(
  type = "Animal",
  yi = as.numeric(animal_mv_adjusted_model$beta),
  ci.lb = as.numeric(animal_mv_adjusted_model$ci.lb),
  ci.ub = as.numeric(animal_mv_adjusted_model$ci.ub),
  label = "Adjusted animal rma.mv model"
)


##################
### Combine original wide data with adjusted estimates for plotting

taar1_animal_final_merge <- taar1_animal %>%
  left_join(
    taar1_animal_final %>% select(result_id, yi_adj, vi_adj),
    by = "result_id"
  )

taar1_human_final_merge <- taar1_human_4 %>%
  left_join(
    taar1_human_final %>% select(result_id, yi_adj, vi_adj),
    by = "result_id"
  )


##################
### Human random-effects summaries for the triangulation plot

### The visible human rows are individual RCTs. The original and adjusted
### subgroup summaries are estimated with meta::metagen() from the study-level
### SMDs and standard errors stored in human_taar1.csv. This keeps the
### triangulation plot aligned with the human meta-analysis input file and
### avoids rebuilding the human dataframe manually in this script.

human_original_meta <- meta::metagen(
  TE = yi,
  seTE = sqrt(vi),
  studlab = study_name_drug,
  data = taar1_human_final_merge,
  sm = "SMD",
  common = TRUE,
  random = TRUE,
  prediction = FALSE,
  method.tau = "REML",
  method.random.ci = "classic"
)

human_adjusted_meta <- meta::metagen(
  TE = yi_adj,
  seTE = sqrt(vi_adj),
  studlab = study_name_drug,
  data = taar1_human_final_merge,
  sm = "SMD",
  common = TRUE,
  random = TRUE,
  prediction = FALSE,
  method.tau = "REML",
  method.random.ci = "classic"
)

human_external_original_summary <- tibble::tibble(
  type = "RCT",
  yi = as.numeric(human_original_meta$TE.random),
  ci.lb = as.numeric(human_original_meta$lower.random),
  ci.ub = as.numeric(human_original_meta$upper.random),
  label = "Human metagen random-effects model (4 RCTs)"
)

human_external_adjusted_summary <- tibble::tibble(
  type = "RCT",
  yi = as.numeric(human_adjusted_meta$TE.random),
  ci.lb = as.numeric(human_adjusted_meta$lower.random),
  ci.ub = as.numeric(human_adjusted_meta$upper.random),
  label = "Adjusted human model"
)

external_original_summaries <- bind_rows(
  animal_external_original_summary,
  human_external_original_summary
)

external_adjusted_summaries <- bind_rows(
  animal_external_adjusted_summary,
  human_external_adjusted_summary
)

taar1_plot_all <- bind_rows(
  taar1_animal_final_merge,
  taar1_human_final_merge
)


##################
### Meta-analysis models for explicit inspection

model_unadjusted_all <- metafor::rma(
  yi = yi,
  vi = vi,
  data = taar1_plot_all,
  slab = taar1_plot_all$study,
  method = "DL"
)

model_adjusted_all <- metafor::rma(
  yi = yi_adj,
  vi = vi_adj,
  data = taar1_plot_all,
  slab = taar1_plot_all$study,
  method = "DL"
)


##################
### Plot animal and human studies together

### rob_direction() is sourced above from Functions/BiasCorrect_robviz.R. It is
### an amended copy of the package plotting internals that displays D1-D10. The
### package's exported tri_plot_bias_direction() still draws only D1-D7.
pdf("TAAR1_triangulate_plot.pdf", width = 18, height = 9)
rob_direction(
  taar1_plot_all,
  grouping_levels = c("Animal", "RCT"),
  title = "TAAR1 agonists: animal and human evidence",
  atransf = identity,
  external_subgroup_summaries = external_original_summaries,
  external_subgroup_summaries_adj = external_adjusted_summaries
)
dev.off()
