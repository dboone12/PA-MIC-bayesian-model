install.packages("brms")
install.packages("rstan")
library(rstan)
stan_version()
options(brms.backend = "rstan")
library(brms)
library(ggplot2)
library(readxl)
mic <- read_excel("~/Desktop/Work Stuff/MIC/excel long format data/Pa_MIC_combined.xlsx")

#prepare variables
mic$strain <- factor(mic$strain)
mic$operator <- factor(mic$operator)
mic$run_id <- factor(mic$run_id)
mic$replicate <- factor(mic$replicate)
mic$antibiotic <- factor(mic$antibiotic)
mic$mic_ug_ml <- as.numeric(mic$mic_ug_ml)
mic$log2_mic <- log2(mic$mic_ug_ml)
mic$mic_censored <- as.character(mic$mic_censored)

#check which ABX is reference
levels(mic$antibiotic)

#build priors
priors <- c(brms::prior(normal(1.5, 0.75), class = b),
            brms::prior(normal(6, 2), class = Intercept),
            brms::prior(normal(0, 1), class = sd),
            brms::prior(normal(0, 1), class = sigma))
#expecting most of the MICs to fall between 32 - 256 ug/mL, also expecting ETP MIC to be 1 - 2 dilution higher than CAZ

#fit model
set.seed(123)
fit <- brm(
  log2_mic | cens(mic_censored) ~ antibiotic +
    (antibiotic | strain) +
    (1 | operator), data = mic, 
  family = gaussian(),
  prior = priors,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = 0.95),
  seed = 42)
#error
