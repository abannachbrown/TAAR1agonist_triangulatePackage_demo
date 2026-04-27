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
