# TAAR1 agonists - bias-corrected Bayesian meta-analysis using jarbes
#
# This script applies the jarbes bias-correction approach (Verde 2021,
# Biometrical Journal) to the GALENOS TAAR1 agonist datasets. The animal studies
# play the role of "observational" (potentially biased) evidence, and the
# RCTs play the role of the unbiased reference arm, mirroring the setup in
# the ppvipd example shipped with the jarbes package.
#
# Reference: Verde PE (2021).
#   A bias-corrected meta-analysis model for combining, studies of different types and quality.
#   Biometrical Journal 63(2):406-422.
#   https://doi.org/10.1002/bimj.201900376
#
# See also: sandbox/testing_jarbes.R for a worked example with the
# jarbes built-in ppvipd dataset.


# ===========================================================================
# 1. Packages
# ===========================================================================

# Install if needed (comment out once installed):
# install.packages("jarbes")
# install.packages("rjags")
# install.packages("R2jags")
# install.packages("LearnBayes")

library(jarbes)
library(rjags)
library(R2jags)
library(LearnBayes)
library(dplyr)
library(readr)
library(stringr)

# ===========================================================================
# 2. Load animal data (drug-level estimates - multi-level model, n = 17 rows)
# ===========================================================================

# taar1_drug_merged.csv contains one row per drug-level animal estimate with
# SYRCLE-style risk-of-bias judgements.
# Columns used here:
#   DrugName  - drug identifier
#   yi        - SMD point estimate
#   se        - standard error (jarbes expects seTE, not variance)

taar1_animal_raw <- read_csv("taar1_drug_merged.csv")

# Rename to the variable names expected throughout this script.
taar1_animal <- taar1_animal_raw %>%
  rename(name = DrugName) %>%
  mutate(
    design = "Animal",
    TE    = yi,
    seTE  = sqrt(vi)
  ) %>%
  select(name, design, TE, seTE, everything())


# ===========================================================================
# 3. Load human RCT data (trial-level SMDs, n = 4 rows)
# ===========================================================================

# human_taar1.csv contains study-level SMDs (yi), standard errors (se),
# and ROB2 judgements. The jarbes BC model treats these as the unbiased arm.

taar1_human_raw <- read_csv("human_taar1.csv")

taar1_human <- taar1_human_raw %>%
  rename(name = study_name_drug) %>%
  mutate(
    design = "RCT",
    TE    = yi,
    seTE  = sqrt(vi)
  ) %>%
  select(name, design, TE, seTE, everything())


# ===========================================================================
# 4. Combine into a single dataset mirroring the ppvipd structure
# ===========================================================================

# jarbes / the Verde 2021 BUGS script expects a flat data frame with columns:
#   name   - study label
#   design - "RCT" or "Animal" (equivalent to "RCT" vs observational)
#   TE     - treatment effect
#   seTE   - standard error of TE

taar1_combined <- bind_rows(taar1_animal, taar1_human)

# Quick sanity check
cat("Rows by design type:\n")
print(table(taar1_combined$design))
cat("\nFirst few rows:\n")
print(taar1_combined[, c("name", "design", "TE", "seTE")])


# ===========================================================================
# 5. BUGS model: Bayesian random-effects meta-analysis
# ===========================================================================

# This is the standard Bayesian RE model from Verde 2021 (supplementary).
# It is shared by the RCT-only and animal-only analyses in Section 6.
#
# Likelihood:
#   y[i] | theta[i]  ~  Normal(theta[i], se.y[i]^2)   [within-study error known]
#   theta[i] | mu, tau  ~  Normal(mu, tau^2)            [random study effects]
#
# Priors:
#   mu                 ~  Normal(0, 1000)               [diffuse prior on pooled effect]
#   1/tau^2 (prec.tau) ~  scaled-gamma(scale, df)       [Half-Cauchy on tau via dscaled.gamma]
#
# The half-Cauchy scale parameter (scale.sigma.between = 0.5) follows the
# Verde 2021 ppvipd default and is a weakly informative prior for tau on the
# SMD scale. For the animal studies, where observed heterogeneity is larger,
# this prior is still conservative and should not dominate the posterior.
#
# theta.new is the posterior predictive effect for a new study (shrinkage
# estimator), useful for assessing the likely effect of a future study.

cat(
  "model
{
  for (i in 1:N) {
    # Likelihood: observed effect y[i] with known sampling variance se.y[i]^2
    y[i]      ~ dnorm(theta[i], pre.y[i])
    pre.y[i] <- pow(se.y[i], -2)

    # Study-specific random effect drawn from the population distribution
    theta[i]  ~ dnorm(mu, prec.tau)
  }

  # Between-study heterogeneity
  # dscaled.gamma(scale, df) gives a half-Cauchy(scale) marginal prior on tau
  tau        <- 1 / sqrt(prec.tau)
  prec.tau   ~ dscaled.gamma(scale.sigma.between, df.scale.between)

  # Pooled effect (population mean)
  mu         ~ dnorm(0.0, 0.001)

  # Posterior predictive effect for a hypothetical new study
  theta.new  ~ dnorm(mu, prec.tau)
}",
  file = "meta_bayes.bugs"
)

# Shared MCMC settings (Verde 2021 defaults)
n.chains   <- 4
n.iter     <- 50000
n.thin     <- 1
n.burnin   <- 20000

# Half-Cauchy scale for tau: 0.5 SMD units (Verde 2021 default)
scale.sigma.between <- 0.5
df.scale.between    <- 1

par.re <- c("theta", "mu", "theta.new", "tau")

set.seed(2026)

# ===========================================================================
# 6. Bayesian RE meta-analysis: RCTs and animal studies separately
# ===========================================================================

# ---------------------------------------------------------------------------
# 6a. RCT-only model  (n = 4 trials)
# ---------------------------------------------------------------------------

# The RCT estimates are on the SMD scale; all four trials compare a TAAR1
# agonist (ulotaront or ralmitaront) against placebo in adults with acute
# schizophrenia.

y.rct    <- taar1_human$TE
se.y.rct <- taar1_human$seTE
N.rct    <- length(y.rct)

data.rct <- list(
  y                    = y.rct,
  se.y                 = se.y.rct,
  N                    = N.rct,
  scale.sigma.between  = scale.sigma.between,
  df.scale.between     = df.scale.between
)

m.rct <- R2jags::jags(
  data               = data.rct,
  inits              = NULL,
  parameters.to.save = par.re,
  model.file         = "meta_bayes.bugs",
  n.chains           = n.chains,
  n.iter             = n.iter,
  n.thin             = n.thin,
  n.burnin           = n.burnin,
  DIC                = TRUE
)

print(m.rct)

# Extract posterior samples for mu (pooled RCT effect)
attach.jags(m.rct, overwrite = TRUE)
mu.rct <- mu

cat("\n--- RCT pooled SMD ---\n")
cat(sprintf("  Posterior mean:   %.3f\n", mean(mu.rct)))
cat(sprintf("  95%% CrI:         [%.3f, %.3f]\n",
            quantile(mu.rct, 0.025), quantile(mu.rct, 0.975)))


# ---------------------------------------------------------------------------
# 6b. Animal-only model  (n = 17 drug-level estimates)
# ---------------------------------------------------------------------------

# The animal estimates are also on the SMD scale, derived from a multivariate
# meta-analysis of locomotor activity outcomes (see taar1-agonists.R).
# Animal studies play the role of "observational" evidence in the Verde 2021
# framework: informative but potentially systematically biased.

y.animal    <- taar1_animal$TE
se.y.animal <- taar1_animal$seTE
N.animal    <- length(y.animal)

data.animal <- list(
  y                    = y.animal,
  se.y                 = se.y.animal,
  N                    = N.animal,
  scale.sigma.between  = scale.sigma.between,
  df.scale.between     = df.scale.between
)

m.animal <- R2jags::jags(
  data               = data.animal,
  inits              = NULL,
  parameters.to.save = par.re,
  model.file         = "meta_bayes.bugs",
  n.chains           = n.chains,
  n.iter             = n.iter,
  n.thin             = n.thin,
  n.burnin           = n.burnin,
  DIC                = TRUE
)

print(m.animal)

# Extract posterior samples for mu (pooled animal effect)
attach.jags(m.animal, overwrite = TRUE)
mu.animal <- mu

cat("\n--- Animal pooled SMD ---\n")
cat(sprintf("  Posterior mean:   %.3f\n", mean(mu.animal)))
cat(sprintf("  95%% CrI:         [%.3f, %.3f]\n",
            quantile(mu.animal, 0.025), quantile(mu.animal, 0.975)))


# ---------------------------------------------------------------------------
# 6c. Side-by-side comparison of RCT vs animal posterior distributions
# ---------------------------------------------------------------------------

cat("\n========================================\n")
cat("  Comparison: RCT vs Animal RE models\n")
cat("========================================\n")

rct_mean    <- mean(mu.rct)
rct_cri     <- quantile(mu.rct,    probs = c(0.025, 0.975))
animal_mean <- mean(mu.animal)
animal_cri  <- quantile(mu.animal, probs = c(0.025, 0.975))

cat(sprintf("  RCT    mu:  %.3f  [%.3f, %.3f]\n",
            rct_mean,    rct_cri[1],    rct_cri[2]))
cat(sprintf("  Animal mu:  %.3f  [%.3f, %.3f]\n",
            animal_mean, animal_cri[1], animal_cri[2]))

# Interval widths (a measure of estimation uncertainty)
cat(sprintf("\n  RCT    95%% CrI width: %.3f\n", diff(rct_cri)))
cat(sprintf("  Animal 95%% CrI width: %.3f\n", diff(animal_cri)))

# Probability that the animal effect exceeds the RCT effect
# (positive values = animal estimates are larger, consistent with known
# translational inflation)
mu.rct_resamp <- sample(mu.rct, length(mu.animal), replace = TRUE)
p_animal_gt_rct <- mean(mu.animal > mu.rct_resamp)
cat(sprintf("\n  P(animal mu > RCT mu): %.3f\n", p_animal_gt_rct))


# ===========================================================================
# 7. BUGS model: Bias-Corrected (BC) meta-analysis
# ===========================================================================

# The BC model (Verde 2021, Section 2.3) extends the RE model by introducing
# a latent indicator T[i] that allocates each study to one of two components:
#
#   Component 1 (T[i] = 1, "unbiased"):
#     theta[i]      ~ Normal(mu[1], tau^2)
#     Represents studies whose effects are exchangeable with the RCTs.
#
#   Component 2 (T[i] = 2, "biased"):
#     theta.bias[i] ~ Normal(mu[2], tau^2 / w[2,i])   [slash distribution]
#     mu[2]          = mu[1] + B   where B > 0
#     w[2,i]        ~ Beta(nu=0.5, 1)   [heavier tails than Component 1]
#     Represents studies whose effects are inflated by a systematic bias B.
#
# The observed effect is the mixture:
#   theta.bc[i] = theta[i] * (1 - I[i])  +  theta.bias[i] * I[i]
#   where I[i] = T[i] - 1  (0 = unbiased, 1 = biased)
#
# Priors:
#   p.bias[2] ~ Beta(alpha.bias, beta.bias)   [fraction of biased studies]
#   B         ~ Uniform(0, B.max)             [bias is positive: animal > RCT]
#   mu[1]     ~ Normal(0, 100)                [diffuse prior on true pooled effect]
#   tau       ~ half-Cauchy(0.5) via dscaled.gamma
#
# B.max is set to 3 SMD units. The observed gap between the animal pooled
# mean (~1.11) and the RCT pooled mean (~0.15) is about 0.96 SMD, so B.max = 3
# is conservative and avoids truncating the prior near the observed difference.
#
# IMPORTANT - sorting of the data:
# Studies are sorted by ascending TE before passing to JAGS. T[1] is fixed to
# component 1 (the smallest TE, likely an RCT) and T[N] to component 2 (the
# largest TE, likely an animal study). This initialisation is required by the
# JAGS dcat sampler to break label switching and identify the two components.

cat(
  "model
{
  for (i in 1:N) {
    # Likelihood: observed effect with known sampling variance
    y[i]          ~ dnorm(theta.bc[i], pre.y[i])
    pre.y[i]      <- pow(se.y[i], -2)

    # Mixture: unbiased component (I=0) or biased component (I=1)
    theta.bc[i]   <- theta[i] * (1 - I[i]) + theta.bias[i] * I[i]
    I[i]          <- T[i] - 1

    # Component-specific random effects
    theta[i]      ~ dnorm(mu[1], prec.tau[i])
    theta.bias[i] ~ dnorm(mu[2], prec.tau[i])

    # Latent study allocation
    T[i]          ~ dcat(p.bias[1:2])

    # Slash parameterisation: heavier tails for biased component
    prec.tau[i]   <- inv.var[T[i]] * w[T[i], i]
    w[1, i]       <- 1
    w[2, i]       ~ dbeta(nu, 1)
  }

  nu <- 1/2   # slash shape parameter (Verde 2021)

  # Prior on fraction of biased studies
  p.bias[2] ~ dbeta(alpha.bias, beta.bias)
  p.bias[1] <- 1 - p.bias[2]

  # Between-study heterogeneity (shared across components)
  tau          <- 1 / sqrt(inv.var[1])
  inv.var[1]   ~ dscaled.gamma(scale.sigma.between, df.scale.between)
  inv.var[2]   <- inv.var[1]

  # Pooled effects
  mu[1] ~ dnorm(0.0, 0.01)      # bias-corrected (unbiased component) mean
  B     ~ dunif(0, B.max)       # systematic bias: B > 0 means animal > RCT
  mu[2] <- mu[1] + B            # biased component mean
}",
  file = "BC.bugs"
)


# ===========================================================================
# 8. Fit the BC model to the combined TAAR1 dataset
# ===========================================================================

# ---------------------------------------------------------------------------
# 8a. Prepare combined data (sorted by TE, as required by the BC model)
# ---------------------------------------------------------------------------

# Sort all 21 rows (17 animal + 4 RCT) by ascending effect estimate.
# After sorting, T[1] is fixed to component 1 (smallest TE, unbiased anchor)
# and T[N] to component 2 (largest TE, biased anchor).

taar1_sorted <- taar1_combined %>%
  arrange(TE)

y.bc    <- taar1_sorted$TE
se.y.bc <- taar1_sorted$seTE
N.bc    <- nrow(taar1_sorted)

# T vector: NA for free studies; T[1]=1 anchors the smallest effect to the
# unbiased component, T[N]=2 anchors the largest to the biased component.
T.init      <- rep(NA, N.bc)
T.init[1]   <- 1
T.init[N.bc] <- 2

# ---------------------------------------------------------------------------
# 8b. Prior on p.bias: Beta(alpha.bias, beta.bias)
# ---------------------------------------------------------------------------

# The prior encodes the expectation that roughly N.obs / N.total studies are
# biased (i.e. the animal studies), with some probability that RCTs also fall
# into the biased component.
#
# We follow Verde 2021: place the prior median at (N.obs - 1) / N and the
# 90th percentile at N.obs / N. This concentrates mass below N.obs / N while
# allowing a small tail above (some RCTs could be biased; some animal studies
# unbiased).

N.ani.bc <- nrow(taar1_animal)   # 17 animal studies
N.bc                             # 21 total

median.prior   <- (N.ani.bc - 1) / N.bc   # 16/21 ≈ 0.762
percentile90   <- N.ani.bc / N.bc         # 17/21 ≈ 0.810

cat(sprintf("\nPrior on p.bias:  median = %.3f,  90th pctile = %.3f\n",
            median.prior, percentile90))

beta.par   <- LearnBayes::beta.select(
  list(p = 0.50, x = median.prior),
  list(p = 0.90, x = percentile90)
)
alpha.bias <- beta.par[1]
beta.bias  <- beta.par[2]

cat(sprintf("Beta prior parameters: alpha = %.3f,  beta = %.3f\n",
            alpha.bias, beta.bias))

# ---------------------------------------------------------------------------
# 8c. Run the BC model
# ---------------------------------------------------------------------------

B.max <- 3   # maximum plausible bias on the SMD scale (see Section 7 note)

data.bc <- list(
  y                   = y.bc,
  se.y                = se.y.bc,
  N                   = N.bc,
  T                   = T.init,
  scale.sigma.between = scale.sigma.between,
  df.scale.between    = df.scale.between,
  alpha.bias          = alpha.bias,
  beta.bias           = beta.bias,
  B.max               = B.max
)

par.bc <- c(
  "theta.bc",     # posterior study-level bias-corrected effects
  "theta",        # unbiased component study effects
  "theta.bias",   # biased component study effects
  "mu",           # mu[1] = bias-corrected pooled SMD; mu[2] = biased mean
  "tau",          # between-study SD
  "p.bias",       # p.bias[2] = posterior fraction of biased studies
  "I",            # I[i] = 1 if study i allocated to biased component
  "w",            # slash weights for biased component
  "B"             # estimated systematic bias (B = mu[2] - mu[1])
)

set.seed(2026)

bc.taar1 <- R2jags::jags(
  data               = data.bc,
  inits              = NULL,
  parameters.to.save = par.bc,
  model.file         = "BC.bugs",
  n.chains           = n.chains,
  n.iter             = 50000,
  n.thin             = 4,       # thin more due to larger model
  n.burnin           = 20000,
  DIC                = TRUE
)

print(bc.taar1)

# ---------------------------------------------------------------------------
# 8d. Extract key posterior quantities
# ---------------------------------------------------------------------------

attach.jags(bc.taar1, overwrite = TRUE)

# mu[,1] = bias-corrected pooled SMD (unbiased component mean)
# mu[,2] = biased component mean (should be close to the animal pooled effect)
mu.bc.corrected <- mu[, 1]
mu.bc.biased    <- mu[, 2]
B.posterior     <- B
p.bias.post     <- p.bias[, 2]


# ---------------------------------------------------------------------------
# 8e. Summary of BC model results
# ---------------------------------------------------------------------------

cat("\n=================================================\n")
cat("  BC model: key posterior summaries\n")
cat("=================================================\n")

cat(sprintf("\n  mu[1]  (bias-corrected pooled SMD):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(mu.bc.corrected)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(mu.bc.corrected, 0.025),
            quantile(mu.bc.corrected, 0.975)))

cat(sprintf("\n  mu[2]  (biased component mean):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(mu.bc.biased)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(mu.bc.biased, 0.025),
            quantile(mu.bc.biased, 0.975)))

cat(sprintf("\n  B  (systematic bias = mu[2] - mu[1]):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(B.posterior)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(B.posterior, 0.025),
            quantile(B.posterior, 0.975)))

cat(sprintf("\n  p.bias  (posterior fraction of biased studies):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(p.bias.post)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(p.bias.post, 0.025),
            quantile(p.bias.post, 0.975)))

# Effective number of unbiased studies
eff.unbiased <- mean((1 - p.bias.post) * N.bc)
cat(sprintf("\n  Effective number of unbiased studies: %.1f  (of %d total)\n",
            eff.unbiased, N.bc))


# ---------------------------------------------------------------------------
# 8f. Study-level posterior probability of being biased
# ---------------------------------------------------------------------------

# I[i] is the posterior mean of the biased indicator for each sorted study.
# Values near 1 indicate the model strongly classifies a study as biased.

I.post <- apply(I, 2, mean)   # I is [n.samples x N.bc]

study_bias_summary <- taar1_sorted %>%
  select(name, design, TE, seTE) %>%
  mutate(
    P_biased = round(I.post, 3)
  ) %>%
  arrange(desc(P_biased))

cat("\n  Study-level P(biased) - top 10:\n")
print(head(study_bias_summary, 10))


# ===========================================================================
# 9. Three-way comparison: RCT-only, Animal-only, BC model
# ===========================================================================

cat("\n=================================================\n")
cat("  Final comparison across models\n")
cat("=================================================\n")
cat("  (all values are posterior mean [95% CrI])\n\n")

fmt <- function(samps) {
  sprintf("%.3f  [%.3f, %.3f]",
          mean(samps),
          quantile(samps, 0.025),
          quantile(samps, 0.975))
}

cat(sprintf("  RCT-only mu:             %s\n", fmt(mu.rct)))
cat(sprintf("  Animal-only mu:          %s\n", fmt(mu.animal)))
cat(sprintf("  BC mu[1] (corrected):    %s\n", fmt(mu.bc.corrected)))
cat(sprintf("  BC B (bias):             %s\n", fmt(B.posterior)))

# Width of 95% CrI for mu[1] vs RCT-only (informativeness gain from borrowing)
width_rct <- diff(quantile(mu.rct,          probs = c(0.025, 0.975)))
width_bc  <- diff(quantile(mu.bc.corrected, probs = c(0.025, 0.975)))

cat(sprintf("\n  95%% CrI width: RCT-only = %.3f,  BC mu[1] = %.3f\n",
            width_rct, width_bc))
cat(sprintf("  Precision gain from BC model: %.1fx narrower\n",
            width_rct / width_bc))

# Save named copies before attach.jags() overwrites mu etc. again
mu.bc.corrected.inf <- mu.bc.corrected
B.posterior.inf     <- B.posterior
p.bias.post.inf     <- p.bias.post


# ===========================================================================
# 10. BC model with non-informative prior on p.bias
#     ("data-driven" version, following Verde 2021 Section 3.2.2)
# ===========================================================================

# In the stemcells example Verde uses alpha.bias = beta.bias = 1 (a uniform
# prior on the fraction of biased studies). This removes any assumption about
# how many studies are biased and lets the data alone identify the mixture.
# For TAAR1 this answers: can the observed effect sizes alone separate the
# two components, without us telling the model that 17/21 are likely biased?
#
# The same BC.bugs model is reused; only the hyperparameters change.

data.bc.noninf <- data.bc   # copy sorted y, se.y, N, T, tau priors, B.max
data.bc.noninf$alpha.bias <- 1
data.bc.noninf$beta.bias  <- 1

set.seed(2026)

bc.taar1.noninf <- R2jags::jags(
  data               = data.bc.noninf,
  inits              = NULL,
  parameters.to.save = par.bc,
  model.file         = "BC.bugs",
  n.chains           = n.chains,
  n.iter             = 50000,
  n.thin             = 4,
  n.burnin           = 20000,
  DIC                = TRUE
)

print(bc.taar1.noninf)

attach.jags(bc.taar1.noninf, overwrite = TRUE)

mu.bc.noninf    <- mu[, 1]
B.noninf        <- B
p.bias.noninf   <- p.bias[, 2]

cat("\n=================================================\n")
cat("  BC model (non-informative prior): summaries\n")
cat("=================================================\n")

cat(sprintf("\n  mu[1]  (bias-corrected pooled SMD):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(mu.bc.noninf)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(mu.bc.noninf, 0.025),
            quantile(mu.bc.noninf, 0.975)))

cat(sprintf("\n  B  (systematic bias):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(B.noninf)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(B.noninf, 0.025),
            quantile(B.noninf, 0.975)))

cat(sprintf("\n  p.bias  (posterior fraction of biased studies):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(p.bias.noninf)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(p.bias.noninf, 0.025),
            quantile(p.bias.noninf, 0.975)))

# Probability that a substantial bias component exists in the data
# (Verde 2021 diagnostic: P(B > 2*tau))
tau.noninf    <- tau
delta.noninf  <- mu[, 2] - mu[, 1]   # = B posterior samples
cut.point.ni  <- 2 * mean(tau.noninf)
p.bias.exists <- mean(delta.noninf > cut.point.ni)
cat(sprintf("\n  P(B > 2*tau) - probability a bias component exists: %.3f\n",
            p.bias.exists))

# Save named copies
mu.bc.noninf.saved  <- mu.bc.noninf
B.noninf.saved      <- B.noninf
p.bias.noninf.saved <- p.bias.noninf


# ===========================================================================
# 11. PBias model: per-study probability of bias linked to ROB score
#     (Verde 2021 Section 3.2.6, heart disease / bone marrow example)
# ===========================================================================

# ---------------------------------------------------------------------------
# 11a. Derive a scalar ROB covariate x[i] for each study
# ---------------------------------------------------------------------------

# The stemcells example used "number of discrepancies" (Cochrane ROB1 tool)
# as the covariate x[i]. For TAAR1 we derive an analogous score from the
# available domain judgements, computed separately per evidence type because
# the two arms use different column names and different value codings.
#
# Animal studies (SYRCLE, taar1_drug_merged.csv):
#   Columns: "D1 Allocation sequence" ... "D10 Free from other risks of bias"
#   Values:  "Yes" = low risk (good); "No" = high risk (bad);
#            "Unclear" = moderate (treated as non-low)
#   Score:   proportion of the 10 domains coded "No" or "Unclear"
#
# Human RCTs (ROB2, human_taar1.csv):
#   Columns: "Bias due to randomization" ... "Bias due to period and carryover effects"
#   Values:  "Low" = good; "Some concerns" = moderate; "High" = bad;
#            "NA" (string) = not applicable (excluded from denominator)
#   Score:   proportion of applicable domains coded "Some concerns" or "High"
#
# Both scores are on the 0-1 scale, so the logistic slope alpha.1 is
# interpretable uniformly across arms.

# -- Animal ROB score -------------------------------------------------------
animal_rob_cols <- c(
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
)

animal_rob_scores <- taar1_animal %>%
  select(name, all_of(animal_rob_cols)) %>%
  rowwise() %>%
  mutate(
    rob_score = {
      vals <- c_across(all_of(animal_rob_cols))
      # "No" = high risk (1), "Unclear" = half-point (0.5), "Yes" = low risk (0)
      mean(dplyr::case_when(
        vals == "No"      ~ 1,
        vals == "Unclear" ~ 0.5,
        vals == "Yes"     ~ 0,
        TRUE              ~ NA_real_
      ), na.rm = TRUE)
    }
  ) %>%
  ungroup() %>%
  select(name, rob_score)

# -- Human ROB score --------------------------------------------------------
human_rob_cols <- c(
  "Bias due to randomization",
  "Bias due to deviations from intended intervention",
  "Bias due to missing data",
  "Bias due to outcome measurement",
  "Bias due to selected reported results",
  "Bias due to period and carryover effects"
)

human_rob_scores <- taar1_human %>%
  select(name, all_of(human_rob_cols)) %>%
  rowwise() %>%
  mutate(
    rob_score = {
      vals       <- c_across(all_of(human_rob_cols))
      applicable <- vals[!is.na(vals) & vals != "NA"]
      if (length(applicable) == 0) NA_real_
      else mean(applicable %in% c("Some concerns", "High"))
    }
  ) %>%
  ungroup() %>%
  select(name, rob_score)

# -- Join scores into taar1_sorted ------------------------------------------
rob_lookup <- bind_rows(animal_rob_scores, human_rob_scores)

taar1_sorted <- taar1_sorted %>%
  select(-any_of("rob_score")) %>%        # drop the previous all-NA attempt
  left_join(rob_lookup, by = "name")

cat("\nROB scores (proportion of non-Low / non-Yes domains) by study:\n")
print(
  taar1_sorted %>%
    select(name, design, TE, rob_score) %>%
    arrange(design, rob_score),
  n = Inf
)

# x[i] vector in the same sorted-by-TE order used by y.bc
x.rob <- taar1_sorted$rob_score


# ---------------------------------------------------------------------------
# 11b. BUGS model: PBias - per-study bias probability via logistic regression
# ---------------------------------------------------------------------------

# Instead of a single global p.bias[2] ~ Beta(...), the PBias model gives
# each study its own probability of being in the biased component:
#
#   logit(P[i, 2]) = alpha.0 + alpha.1 * x[i]
#
# where x[i] is the ROB score for study i.
#
# alpha.0: baseline log-odds of bias when ROB score = 0 (all domains Low)
# alpha.1: change in log-odds per unit increase in ROB score
#          positive alpha.1 => higher ROB score => higher probability of bias
#
# Priors (Verde 2021 defaults, weakly informative on logit scale):
#   alpha.0 ~ Normal(0, precision=0.1)   [SD ~ 3.16]
#   alpha.1 ~ Normal(0, precision=0.1)   [SD ~ 3.16]
#
# With a 0-1 ROB covariate, a prior SD of ~3 on alpha.1 allows the
# log-odds to swing by ±3 across the full covariate range -- large enough
# for the posterior to be data-driven but not completely uninformative.

cat(
  "model
{
  for (i in 1:N) {
    # Likelihood
    y[i]        ~ dnorm(theta.bc[i], pre.y[i])
    pre.y[i]   <- pow(se.y[i], -2)

    # Mixture: unbiased (T=1) or biased (T=2) component
    theta.bc[i] <- theta[i] * (1 - I[i]) + theta.bias[i] * I[i]
    I[i]        <- T[i] - 1

    theta[i]      ~ dnorm(mu[1], prec.tau[i])
    theta.bias[i] ~ dnorm(mu[2], prec.tau[i])

    # Study-specific allocation probability from ROB score
    T[i] ~ dcat(P[i, 1:2])

    logit(P[i, 2]) <- alpha.0 + alpha.1 * x[i]
    P[i, 1]        <- 1 - P[i, 2]

    # Slash parameterisation for biased component
    prec.tau[i] <- inv.var[T[i]] * w[T[i], i]
    w[1, i]     <- 1
    w[2, i]     ~ dbeta(nu, 1)
  }

  nu <- 1/2

  # Logistic regression priors (weakly informative on logit scale)
  alpha.0 ~ dnorm(prior.mean.alpha.0, prior.scale.alpha.0)
  alpha.1 ~ dnorm(prior.mean.alpha.1, prior.scale.alpha.1)

  # Between-study heterogeneity
  tau        <- 1 / sqrt(inv.var[1])
  inv.var[1] ~ dscaled.gamma(scale.sigma.between, df.scale.between)
  inv.var[2] <- inv.var[1]

  # Pooled effects
  mu[1] ~ dnorm(0.0, 0.01)
  B     ~ dunif(0, B.max)
  mu[2] <- mu[1] + B
}",
  file = "PBias.bugs"
)


# ---------------------------------------------------------------------------
# 11c. Fit the PBias model
# ---------------------------------------------------------------------------

data.pbias <- list(
  y                    = y.bc,
  se.y                 = se.y.bc,
  x                    = x.rob,
  N                    = N.bc,
  T                    = T.init,
  scale.sigma.between  = scale.sigma.between,
  df.scale.between     = df.scale.between,
  B.max                = B.max,
  prior.mean.alpha.0   = 0,
  prior.scale.alpha.0  = 0.1,   # precision = 0.1  =>  SD ~ 3.16
  prior.mean.alpha.1   = 0,
  prior.scale.alpha.1  = 0.1
)

par.pbias <- c(
  "theta.bc",   # study-level bias-corrected effects
  "mu",         # mu[1] = bias-corrected pooled SMD
  "tau",        # between-study SD
  "P",          # P[i, 2] = per-study probability of being biased
  "T",          # latent study allocation
  "w",          # slash weights
  "B",          # systematic bias
  "alpha.0",    # logistic intercept
  "alpha.1"     # logistic slope on ROB score
)

set.seed(2026)

PBias.taar1 <- R2jags::jags(
  data               = data.pbias,
  inits              = NULL,
  parameters.to.save = par.pbias,
  model.file         = "PBias.bugs",
  n.chains           = n.chains,
  n.iter             = 50000,
  n.thin             = 4,
  n.burnin           = 20000,
  DIC                = TRUE
)

print(PBias.taar1)


# ---------------------------------------------------------------------------
# 11d. Extract PBias posteriors
# ---------------------------------------------------------------------------

attach.jags(PBias.taar1, overwrite = TRUE)

mu.pbias      <- mu[, 1]
B.pbias       <- B
alpha.0.post  <- alpha.0
alpha.1.post  <- alpha.1

cat("\n=================================================\n")
cat("  PBias model: key posterior summaries\n")
cat("=================================================\n")

cat(sprintf("\n  mu[1]  (bias-corrected pooled SMD):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(mu.pbias)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(mu.pbias, 0.025),
            quantile(mu.pbias, 0.975)))

cat(sprintf("\n  B  (systematic bias):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(B.pbias)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(B.pbias, 0.025),
            quantile(B.pbias, 0.975)))

cat(sprintf("\n  alpha.0  (logistic intercept):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(alpha.0.post)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(alpha.0.post, 0.025),
            quantile(alpha.0.post, 0.975)))

cat(sprintf("\n  alpha.1  (logistic slope on ROB score):\n"))
cat(sprintf("    Mean:   %.3f\n",   mean(alpha.1.post)))
cat(sprintf("    95%% CrI: [%.3f, %.3f]\n",
            quantile(alpha.1.post, 0.025),
            quantile(alpha.1.post, 0.975)))

# Posterior probability that alpha.1 > 0
# (positive = higher ROB score associated with higher probability of being biased)
p_alpha1_pos <- mean(alpha.1.post > 0)
cat(sprintf("\n  P(alpha.1 > 0): %.3f\n", p_alpha1_pos))
cat("  [positive = higher ROB score predicts higher probability of being biased]\n")


# ---------------------------------------------------------------------------
# 11e. Per-study posterior probability of bias from PBias model
# ---------------------------------------------------------------------------

# P[i, 2] is now study-specific (driven by x[i] and alpha.0/alpha.1).
# Compare with the study-level I[i] from the informative BC model (Section 8).

# P array from JAGS: [n.samples x N x 2]; we want column 2 (biased)
P.post.mean <- apply(P[, , 2], 2, mean)

study_pbias_summary <- taar1_sorted %>%
  select(name, design, TE, rob_score) %>%
  mutate(
    P_biased_PBias = round(P.post.mean, 3)
  ) %>%
  arrange(desc(P_biased_PBias))

cat("\n  Per-study P(biased) from PBias model:\n")
print(study_pbias_summary)


# ---------------------------------------------------------------------------
# 11f. Full four-model comparison
# ---------------------------------------------------------------------------

cat("\n=================================================\n")
cat("  Four-model comparison\n")
cat("=================================================\n")
cat("  (all values: posterior mean  [95% CrI])\n\n")

cat(sprintf("  RCT-only RE:              %s\n", fmt(mu.rct)))
cat(sprintf("  Animal-only RE:           %s\n", fmt(mu.animal)))
cat(sprintf("  BC (informative prior):   %s\n", fmt(mu.bc.corrected.inf)))
cat(sprintf("  BC (non-informative):     %s\n", fmt(mu.bc.noninf.saved)))
cat(sprintf("  PBias (ROB covariate):    %s\n", fmt(mu.pbias)))
