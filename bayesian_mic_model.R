library(brms)
library(readxl)
library(tidymodels)
library(skimr)


# data prep ---------------------------------------------------------------

mic <- read_excel("~/Desktop/Work Stuff/MIC/excel long format data/Pa_MIC_combined.xlsx")

skim(mic)

mic <- mic %>% mutate(
  log2mic = log2(mic_ug_ml),
  censored = ifelse(mic_censored == "right", "right", "none"),
  antibiotic = factor(antibiotic, levels = c("CAZ", "ETP")),
  operator = factor(operator),
  run_id = factor(run_id),
  strain = factor(strain)
)
#log2 transformation of the outcome variable to linearizes the doubling dilution series nature of mic assays
#CAZ is setup as the reference level, model intercepts will represent CAZ grand mean

priors <- c(
  prior(normal(5.6, 0.8), class = Intercept),
  prior(normal(1.5, 0.5), class = b, coef = antibioticETP),
  prior(exponential(1), class = sd, group = strain),
  prior(exponential(1), class = sd, group = run_id),
  prior(exponential(1), class = sd, group = operator),
  prior(exponential(1), class = sigma)
)
#Priors
#Intercept N(5.6, 0.8) I expect grand CAZ mean to be between 32 and 64 ug, 5.6 represents the log2 transformation of 48 ug/ml. 
#ETP expected to be approximately 1 - 2 dilutions higher than CAZ on log2 scale
#standard deviations must be positive, hence exponential distribution

fit <- brm(
  log2mic | cens(censored) ~ antibiotic + (1 | strain) + (1 | run_id) + (1 | operator),
  data = mic,
  family = gaussian(),
  prior = priors,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 42,
  backend = "cmdstanr"
) 

#Model 1 - random intercepts only
#log2mic is the outcome cens(cesored) treats right censored observations as censored (<256ug/ml) rather than exactly 256 ug/ml
#antibiotic is a fixed effect - the effect we want to estimate
#(1 | strain) random effect of strain
#(1 | run_id) random effect of day to day variation
#(1 | operator) random effect of operator variability
             
summary(fit)           

plot(fit)

#install.packages("tidybayes")
library(tidybayes)

fit %>%
  spread_draws(r_strain[strain, term]) %>%
  filter(term == "Intercept") %>%
  median_qi(.width = c(0.95)) %>%
  ggplot(aes(y = reorder(strain, r_strain), x = r_strain, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Strain-level deviation from grand mean (log2 scale)",
       y = "Strain",
       title = "Strain random effects") +
  theme_bw()
#plot of posterior median deviation from grand mean after accounting for random effects

ggsave("strain_random_effects.png", width = 10, height = 8, dpi = 300)

pp_check(fit, ndraws = 100)
#dark line is our oberved data
#light blue lines are 100 simulated from fitted model
#our data appears bimodal since CAZ has a lower mic on average than ETP, the fitted model "smooths" this out because both antibiotics are treated as single gaussian family.


library(dplyr)

fit %>%
  spread_draws(r_strain[strain, term], b_Intercept) %>%
  filter(term == "Intercept") %>%
  mutate(strain_mean_log2 = b_Intercept + r_strain,
         strain_mean_mic = 2^strain_mean_log2) %>%
  group_by(strain) %>%
  median_qi(strain_mean_mic, .width = 0.95) %>%
  arrange(desc(strain_mean_mic)) %>%
  print(n = 29)
#estimated MICs for CAZ 

# updated model -----------------------------------------------------------
priors2 <- c(
  prior(normal(5.6, 0.8), class = Intercept),
  prior(normal(1.5, 0.5), class = b, coef = antibioticETP),
  prior(exponential(1), class = sd, group = strain),
  prior(exponential(1), class = sd, group = run_id),
  prior(exponential(1), class = sd, group = operator),
  prior(exponential(1), class = sigma),
  prior(lkj(2), class = cor, group = strain)
)

fit2 <- brm(
  log2mic | cens(censored) ~ antibiotic + (antibiotic | strain) + (1 | run_id) + (1 | operator),
  data = mic,
  family = gaussian(),
  prior = priors2,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 42,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95)
)
#Model 2 - random slope for antibiotic
#family = Gaussian not changed - bimodal distribution not addressed yet?

summary(fit2)

fit2 %>%
  spread_draws(r_strain[strain, term]) %>%
  filter(term == "antibioticETP") %>%
  group_by(strain) %>%
  median_qi(r_strain, .width = 0.95) %>%
  arrange(desc(r_strain)) %>%
  ggplot(aes(y = reorder(strain, r_strain), x = r_strain, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Strain-specific ETP vs CAZ difference (log2 scale)",
       y = "Strain",
       title = "How much higher is ETP MIC vs CAZ MIC per strain?") +
  theme_bw()


fit2 %>%
  spread_draws(r_strain[strain, term], b_Intercept, b_antibioticETP) %>%
  pivot_wider(names_from = term, values_from = r_strain) %>%
  mutate(
    caz_log2 = b_Intercept + Intercept,
    etp_log2 = b_Intercept + b_antibioticETP + Intercept + antibioticETP,
    caz_mic = 2^caz_log2,
    etp_mic = 2^etp_log2
  ) %>%
  group_by(strain) %>%
  median_qi(caz_mic, etp_mic, .width = 0.95) %>%
  arrange(desc(caz_mic)) %>%
  print(n = 29)
#median posteriors for each strain, for each of the two antibiotics along with CI
#back transformed into ug/ml scale

#Model 3 - addresses bimodal distribution
priors3 <- c(
  prior(normal(5.6, 0.8), class = Intercept),
  prior(normal(1.5, 0.5), class = b, coef = antibioticETP),
  prior(exponential(1), class = sd, group = strain),
  prior(exponential(1), class = sd, group = run_id),
  prior(exponential(1), class = sd, group = operator),
  #prior(exponential(1), class = sigma),
  prior(lkj(2), class = cor, group = strain),
  prior(normal(0, 1), class = b, dpar = "sigma")
)

fit3 <- brm(
  bf(log2mic | cens(censored) ~ antibiotic + (antibiotic | strain) + (1 | run_id) + (1 | operator),
     sigma ~ antibiotic),
  data = mic,
  family = gaussian(),
  prior = priors3,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 42,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95)
)

#random effect plot
fit3 %>%
  spread_draws(r_strain[strain, term]) %>%
  filter(term == "Intercept") %>%
  median_qi(r_strain, .width = 0.95) %>%
  ggplot(aes(y = reorder(strain, r_strain), x = r_strain, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Strain-level deviation from grand mean (log2 scale)",
       y = "Strain",
       title = "Strain random effects - fit3") +
  theme_bw()

#posterior medians for all 29 strains
fit3 %>%
  spread_draws(r_strain[strain, term], b_Intercept, b_antibioticETP) %>%
  pivot_wider(names_from = term, values_from = r_strain) %>%
  mutate(
    caz_log2 = b_Intercept + Intercept,
    etp_log2 = b_Intercept + b_antibioticETP + Intercept + antibioticETP,
    caz_mic = 2^caz_log2,
    etp_mic = 2^etp_log2
  ) %>%
  group_by(strain) %>%
  median_qi(caz_mic, etp_mic, .width = 0.95) %>%
  arrange(desc(caz_mic)) %>%
  print(n = 29)


#posterior predictive check
pp_check(fit3, ndraws = 100)

summary(fit3)

plot(fit3)

res <- residuals(fit3)
res_df <- as.data.frame(res)

ggplot(res_df, aes(x = Estimate)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(x = "Residuals", title = "Residual distribution - fit3") +
  theme_bw()
