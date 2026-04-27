# install packages

install.packages("jarbes")
install.packages("rjags")
install.packages("LearnBayes")
library("rjags")
library("jarbes")
library(LearnBayes)


### attempting to follow code presented in Verde 2021 supplementary
### https://onlinelibrary.wiley.com/doi/10.1002/bimj.201900376#bimj2185-bib-0013
### but with Taar-1 data

data("ppvipd")

dat.ipv.all = ppvipd

dat.ipv.all$Risk = c("Unclear", "Low", "Low", "High", "Low", "Low", "High",
                     "Low", "Low","Low", "Low")
dat.ipv.all[, c("name", "design","TE", "seTE", "n.c", "n.v","Risk")]

# set the seed of the random numbers at 2020 (the year of Verde publications)
set.seed(2020)


# Invlogit
invlogit = function(x)
{
  1/(1 + exp(-x))
}


# Bayesian random-effects meta-analysis 

library(R2jags)

# Bayesian Random Effects Model ......................
# Model in BUGS
cat(
  "model
{
for( i in 1 : N ) {
# Likelihood of theta[i] .........................
y[i] ~ dnorm(theta[i], pre.y[i])
pre.y[i] <- pow(se.y[i], -2)
# Simple random effects ..........................
theta[i] ~ dnorm(mu, prec.tau)
}
# Priors for hyper-parameters .......................................
# Variance components ...
tau <- 1/sqrt(prec.tau)
prec.tau ~ dscaled.gamma(scale.sigma.between,
df.scale.between) # spike: tau is small (unbias mu)
# Prior for mu
mu ~ dnorm(0.0, 0.001)
# Posterior predictive distribution
theta.new ~ dnorm(mu, prec.tau)
}",
  file = "meta_bayes.bugs")
###########################################################


# Data model
y = dat.ipv.all$TE
se.y = dat.ipv.all$seTE
N = length(y)
scale.sigma.between = 0.5
df.scale.between = 1

# This list describes the data used by the BUGS script.
data.meta.bayes <- list ("y", "se.y", "N",
                         "scale.sigma.between",
                         "df.scale.between")
# List of parameters
par.meta.bayes <- c("theta",
                    "mu",
                    "theta.new",
                    "tau")


########################################

#### RCTs

y = dat.ipv.all$TE[dat.ipv.all$design=="RCT"]
se.y = dat.ipv.all$seTE[dat.ipv.all$design=="RCT"]
N = length(y)
N

# analyse RCTs only
# Call with OpenBUGS
m.rct <- jags(data.meta.bayes,
              inits = NULL,
              parameters.to.save = par.meta.bayes,
              model = "meta_bayes.bugs",
              n.chains = 4,
              n.iter = 50000,
              n.thin = 1,
              n.burnin = 20000,
              DIC = TRUE)


attach.jags(m.rct, overwrite = TRUE)
mu.rct = mu


########################### 
# Observational Studies

y = dat.ipv.all$TE[dat.ipv.all$design!="RCT"]
se.y = dat.ipv.all$seTE[dat.ipv.all$design!="RCT"]
N = length(y)
N


y.obs = dat.ipv.all$TE[dat.ipv.all$design!="RCT"]
se.y.obs = dat.ipv.all$seTE[dat.ipv.all$design!="RCT"]
N.obs = length(y.obs)
N.obs

# Call with OpenBUGS
m.os <- jags(data.meta.bayes,
             inits = NULL,
             parameters.to.save = par.meta.bayes,
             model = "meta_bayes.bugs",
             n.chains = 4,
             n.iter = 50000,
             n.thin = 1,
             n.burnin = 20000,
             DIC = TRUE)

attach.jags(m.os, overwrite = TRUE)
mu.os = mu


#### Bias Corrected meta-analysis

# The Bias Corrected Model ..........................................
# Model in BUGS
cat(
  "model
{
for( i in 1 : N ) {
# Likelihood of theta[i] ..........................................
y[i] ~ dnorm(theta.bc[i], pre.y[i])
pre.y[i] <- pow(se.y[i], -2)
# Mixture heterocedastic random effects ..........................
theta.bc[i] <- theta[i]*(1-I[i]) + theta.bias[i]*I[i]
I[i] <- T[i] -1
theta[i] ~ dnorm(mu[1], prec.tau[i])
theta.bias[i] ~ dnorm(mu[2], prec.tau[i])
T[i] ~ dcat(p.bias[1:2])
prec.tau[i] <- inv.var[T[i]] * w[T[i],i] #Slash parametrization
w[1,i] <- 1
# slash distribution ...
w[2,i] ~ dbeta(nu, 1)
}
nu <- 1/2
# Priors for hyper-parameters ................................
# Prior of probability classes
p.bias[2] ~ dbeta(alpha.bias, beta.bias)
p.bias[1] <- 1 - p.bias[2]
# Variance components ...
tau <- 1/sqrt(inv.var[1])
inv.var[1] ~ dscaled.gamma(scale.sigma.between,
df.scale.between)
inv.var[2] <- inv.var[1]
# Prior for mu
mu[1] ~ dnorm(0.0, 0.01)
B ~ dunif(0, B.max)
mu[2] <- mu[1] + B
}",
  file = "BC.bugs")

###########################################################
y = sort(dat.ipv.all$TE)
se.y = dat.ipv.all$seTE[order(dat.ipv.all$TE)]
N = length(y)
#T is the index that allocates each study to one of the two components
T = rep(NA, N)
T[1] = 1
T[N] = 2
scale.sigma.between = 0.5
df.scale.between = 1

# probability to be biased
N.obs # number of observational studies 
N # total N

# median of the prior is: (N.obs -1 / N)
median.prior <- (N.obs - 1)/ N # 7/11 or 0.63636
# 90th percentile of the prior is (N.obs / N)
percentile90 <- N.obs / N # 8/11 or 0.727272

# The idea is to limit the number of biased studies up to the number of OS
# and open the possibility to include unbiased studies, which are OS.


library(LearnBayes)


quantile1=list(p = 0.5, x = median.prior)
quantile2=list(p = 0.9, x = percentile90)
beta.par = beta.select(quantile1, quantile2)
alpha.bias = beta.par[1]
beta.bias = beta.par[2]

# Non-informative prior for the probability of bias
# alpha.bias = 1
# beta.bias = 1
# alpha.bias = beta.bias = 1

B.max = 10
# This list describes the data used by the BUGS script.
data.BC <- list ("y", "se.y", "N", "T",
                 "scale.sigma.between",
                 "df.scale.between",
                 "alpha.bias",
                 "beta.bias",
                 "B.max")
# List of parameters
par.BC <- c("theta.bc",
            "theta",
            "theta.bias",
            "mu",
            "tau",
            "p.bias",
            "I",
            "w",
            "B")
# Call with OpenBUGS
bc.1 <- jags(data.BC,
             inits = NULL,
             parameters.to.save = par.BC,
             model = "BC.bugs",
             n.chains = 4,
             n.iter = 50000,
             n.thin = 4,
             n.burnin = 20000,
             DIC = TRUE)

bc.1

attach.jags(bc.1, overwrite = TRUE)
mu.bc = mu[,1]

############################################
install.packages('extraDistr')
library(extraDistr)
# Half-Cauchy distribution 

S = 0.5

# Probability that tau > 0.5 for S = 0.5
1-phcauchy(S, sigma = 0.5)

# Probability that tau > 4*0.5 for S = 0.5
round(1-phcauchy(S*4, sigma = 0.5),3)



########### PPV23 vacine against invasive IPD 

 # Hyper-parameters for πbias ∼ Beta(a0,a1)
alpha.bias
beta.bias

attach.jags(bc.1, overwrite = TRUE)
#### Effective N+ 
N.obs / N


# Posterior probability of bias
attach.jags(bc.1, overwrite = TRUE)
p.bias.1 = p.bias[,2]
round(mean(p.bias.1),2)


# Effective number of studies
round(mean((1-p.bias.1)*N),2)


####### compare Bayesian RE MA and BC Model 

round(quantile(mu.rct, probs = c(0.025, 0.975)),2)

round(mean(mu.rct),2)

RCT.interval = quantile(mu.rct, probs = c(0.025, 0.975))
Length.RCT.interval = RCT.interval[2]- RCT.interval[1]
round(Length.RCT.interval,2)

round(mean(mu.bc),2)

BC.interval = quantile(mu.bc, probs = c(0.025, 0.975))
round(BC.interval, 2)

Length.BC.interval = BC.interval[2]- BC.interval[1]
round(Length.BC.interval,2)


round(Length.RCT.interval/Length.BC.interval*10,2)
