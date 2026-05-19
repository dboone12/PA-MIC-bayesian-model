install.packages("readxl")
library("readxl")
mic <- read_excel("~/Desktop/Work Stuff/Pa_CAZ_MIC.xlsx")
str(mic)
head(mic)

mic$strain <- factor(mic$strain)
mic$operator <- factor(mic$operator)
mic$run_id <- factor(mic$run_id)
mic$replicate <- factor(mic$replicate)
mic$mic_censored <- factor(mic$mic_censored)

mic$mic_ug_ml <- as.numeric(mic$mic_ug_ml)
mic$log2_mic <- log2(mic$mic_ug_ml)

install.packages("ggplot2")
library(ggplot2)

ggplot(mic, aes(x=strain, y=log2_mic, color = operator)) + 
  geom_point(position = position_jitter(width=0.15, height = 0)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust =1 ))

install.packages("lme4")
library(lme4)
fit <- lmer(log2_mic ~ 1 + (1|strain) + (1|operator), data = mic)
summary(fit)

v_strain <- 4.1784
v_op <- 0.1395
v_res <- 1.7117
total <- v_strain + v_op + v_res
v_strain / total
v_op / total 
v_res / total

mu_log2 <- fixef(fit)[1]

re_strain <- ranef(fit)$strain            
strain_ids <- rownames(re_strain)

b_strain <- re_strain[, 1]                

strain_log2_est <- mu_log2 + b_strain
strain_mic_est  <- 2^strain_log2_est

strain_table <- data.frame(
  strain = strain_ids,
  mic_est_ug_ml = as.numeric(strain_mic_est),
  log2_mic_est  = as.numeric(strain_log2_est),
  row.names = NULL
)

strain_table[order(strain_table$mic_est_ug_ml), ]

library(ggplot2)

ggplot(mic, aes(x = strain, y = mic_ug_ml, colour = operator)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.6) +
  geom_point(
    data = strain_table,
    aes(x = strain, y = mic_est_ug_ml),
    inherit.aes = FALSE,
    size = 3,
    shape = 18,
    color = "black"
  ) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Strain",
    y = "MIC (µg/mL, log scale)",
    title = "Raw MIC measurements with model-based strain estimates overlaid"
  )

