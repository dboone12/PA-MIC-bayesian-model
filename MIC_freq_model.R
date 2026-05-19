library(readxl)
library(ggplot2)
library(lme4)
library(tidyr)
mic <- read_excel("~/Desktop/Work Stuff/MIC/excel long format data/Pa_MIC_combined.xlsx")

#prepare variables
mic$strain <- factor(mic$strain)
mic$operator <- factor(mic$operator)
mic$run_id <- factor(mic$run_id)
mic$replicate <- factor(mic$replicate)
mic$antibiotic <- factor(mic$antibiotic)
mic$mic_ug_ml <- as.numeric(mic$mic_ug_ml)
mic$log2_mic <- log2(mic$mic_ug_ml)

#check reference
levels(mic$antibiotic)

ggplot(data = mic, mapping = aes(x = log2_mic)) +
  geom_histogram(bins = 4)

hist(residuals(fit))
qqnorm(residuals(fit))
qqline(residuals(fit))

#fit model
fit <- lmer(log2_mic ~ antibiotic +
              (antibiotic | strain) +
              (1 | operator), 
            data = mic)
summary(fit)

#variance components
var_components <- as.data.frame(VarCorr(fit))
var_components

var_components <- as.data.frame(VarCorr(fit))

v_strain       <- var_components$vcov[1]
v_strain_slope <- var_components$vcov[2]
v_op           <- var_components$vcov[4]
v_res          <- var_components$vcov[5]

total <- v_strain + v_strain_slope + v_op + v_res

cat("Strain intercept variance:   ", round(v_strain / total * 100, 1), "%\n")
cat("Strain slope variance:       ", round(v_strain_slope / total * 100, 1), "%\n")
cat("Operator variance:           ", round(v_op / total * 100, 1), "%\n")
cat("Residual variance:           ", round(v_res / total * 100, 1), "%\n")

#fixed effects
fixef(fit)

#strain level estimates
mu <- fixef(fit)[1]
etp_effect <- fixef(fit)[2]

re_strain <- ranef(fit)$strain
strain_ids <- rownames(re_strain)

#strain estimates for each ABX
strain_CAZ_log2 <- mu + re_strain[, 1]
strain_ETP_log2 <- mu + etp_effect + re_strain[, 1] + re_strain[, 2]

#strain table
strain_table <- data.frame(strain = strain_ids, CAZ_mic_est = 2^strain_CAZ_log2, ETP_mic_est = 2^strain_ETP_log2, CAZ_log2_est = as.numeric(strain_CAZ_log2), ETP_log2_est = as.numeric(strain_ETP_log2), row.names = NULL)

strain_table[order(strain_table$CAZ_mic_est),]


#raw data plot
ggplot(mic, aes(x = strain, y = mic_ug_ml, colour = operator)) +
  geom_point(position=position_jitter(width = 0.15, height = 0), alpha = 0.6) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)) +
  facet_wrap(~ antibiotic, ncol = 1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Strain", 
    y = "MIC (ug/mL, log scale",
    colour = "Operator",
    title = "RAW MIC measurements by strain and antibiotic"
  )

#raw data + strain estimates on top
strain_long <- pivot_longer(strain_table,
                            cols = c(CAZ_mic_est, ETP_mic_est), names_to = "antibiotic", values_to = "mic_est")

strain_long$antibiotic <- ifelse(strain_long$antibiotic == "CAZ_mic_est", "CAZ", "ETP")

ggplot(mic, aes(x = strain, y = mic_ug_ml, colour = operator)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.6) +
  geom_point(
    data = strain_long,
    aes(x = strain, y = mic_est),
    inherit.aes = FALSE,
    size = 3,
    shape = 18,
    color = "black"
  ) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)) +
  facet_wrap(~ antibiotic, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Strain",
    y = "MIC (µg/mL, log scale)",
    colour = "Operator",
    title = "Raw MIC measurements with model-based strain estimates overlaid"
  )

#variance components
var_components <- as.data.frame(VarCorr(fit))
var_components

v_strain   <- var_components$vcov[1]
v_strain_slope <- var_components$vcov[2]
v_op       <- var_components$vcov[4]
v_res      <- var_components$vcov[5]

total <- v_strain + v_strain_slope + v_op + v_res

variance_df <- data.frame(
  group  = c("Strain (Intercept)", "Strain (ETP Slope)", "Operator", "Residual"),
  values = c(v_strain / total, v_strain_slope / total, v_op / total, v_res / total)
)

#bar graph
ggplot(variance_df, aes(x = reorder(group, values), y = values * 100, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, color = "white") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x     = NULL,
    y     = "Percentage of total variance (%)",
    title = "Proportion of variance by source"
  )


#varaince plot lollipop
ggplot(variance_df, aes(x = reorder(group, values), y = values * 100, color = group)) +
  geom_segment(aes(xend = group, y = 0, yend = values * 100), size = 1) +
  geom_point(size = 5) +
  coord_flip() +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x     = NULL,
    y     = "Percentage of total variance (%)",
    title = "Proportion of variance by source"
  )

#ETP vs CAZ estimates
ggplot(strain_table, aes(x = CAZ_mic_est, y = ETP_mic_est, label = strain)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.7, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)) +
  theme_bw() +
  labs(
    x     = "CAZ MIC estimate (µg/mL, log scale)",
    y     = "ETP MIC estimate (µg/mL, log scale)",
    title = "Strain-level CAZ vs ETP MIC estimates"
  )

library(dplyr)
n_distinct(mic$strain)
