# =============================================================================
# MIC Mixed-Effects Model — Ceftazidime (CAZ)
# Pseudomonas aeruginosa | Pa_CAZ_MIC.xlsx
# Model: log2(MIC) ~ 1 + (1|strain) + (1|operator)
# =============================================================================

# --- Load packages ------------------------------------------------------------
if (!requireNamespace("readxl",  quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("lme4",    quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("scales",  quietly = TRUE)) install.packages("scales")

library(readxl)
library(ggplot2)
library(lme4)
library(scales)

# --- Load data ----------------------------------------------------------------
mic <- read_excel("~/Desktop/Work Stuff/MIC/Pa_CAZ_MIC.xlsx")
# str(mic)   # uncomment to inspect structure
# head(mic)  # uncomment to preview data

# --- Factor and numeric conversions ------------------------------------------
mic$strain       <- factor(mic$strain)
mic$operator     <- factor(mic$operator)
mic$run_id       <- factor(mic$run_id)
mic$replicate    <- factor(mic$replicate)
mic$mic_censored <- factor(mic$mic_censored)

mic$mic_ug_ml <- as.numeric(mic$mic_ug_ml)
mic$log2_mic  <- log2(mic$mic_ug_ml)

# Note: mic_censored observations (at detection limit) are currently included.
# If censored values are at a boundary, consider filtering before modelling:
#   mic <- mic[mic$mic_censored == 0, ]

# --- Exploratory plot: raw log2 MIC by strain and operator -------------------
ggplot(mic, aes(x = strain, y = log2_mic, color = operator)) +
  geom_point(position = position_jitter(width = 0.15, height = 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Strain", y = "log2 MIC", title = "Raw CAZ log2 MIC by strain and operator")

# --- Fit mixed-effects model --------------------------------------------------
fit <- lmer(log2_mic ~ 1 + (1 | strain) + (1 | operator), data = mic)
summary(fit)

# --- Extract variance components from model (no manual entry needed) ----------
vc       <- as.data.frame(VarCorr(fit))
v_strain <- vc[vc$grp == "strain",   "vcov"]
v_op     <- vc[vc$grp == "operator", "vcov"]
v_res    <- vc[vc$grp == "Residual", "vcov"]
total    <- v_strain + v_op + v_res

cat(sprintf(
  "Variance components:\n  Strain: %.4f (%.1f%%)\n  Operator: %.4f (%.1f%%)\n  Residual: %.4f (%.1f%%)\n",
  v_strain, 100 * v_strain / total,
  v_op,     100 * v_op     / total,
  v_res,    100 * v_res    / total
))

# --- Compute strain-level MIC estimates --------------------------------------
mu_log2    <- fixef(fit)[1]
re_strain  <- ranef(fit)$strain
strain_ids <- rownames(re_strain)
b_strain   <- re_strain[, 1]

strain_log2_est <- mu_log2 + b_strain
strain_mic_est  <- 2^strain_log2_est

strain_table <- data.frame(
  strain        = strain_ids,
  mic_est_ug_ml = as.numeric(strain_mic_est),
  log2_mic_est  = as.numeric(strain_log2_est),
  row.names     = NULL
)

print(strain_table[order(strain_table$mic_est_ug_ml), ])

# --- Plot: raw data with model estimates overlaid ----------------------------
ggplot(mic, aes(x = strain, y = mic_ug_ml, colour = operator)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.6) +
  geom_point(
    data        = strain_table,
    aes(x = strain, y = mic_est_ug_ml),
    inherit.aes = FALSE,
    size        = 3,
    shape       = 18,
    color       = "black"
  ) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x     = "Strain",
    y     = "MIC (µg/mL, log scale)",
    title = "Raw CAZ MIC measurements with model-based strain estimates overlaid"
  )

# --- Variance components pie chart -------------------------------------------
variance_df <- data.frame(
  group  = c("Strain", "Operator", "Residual"),
  values = c(v_strain / total, v_op / total, v_res / total)
)

ggplot(variance_df, aes(x = "", y = values, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(
    aes(label = paste0(group, "\n", percent(values, accuracy = 0.1))),
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_brewer("Source") +
  labs(title = "Variance Components — CAZ MIC Model")
