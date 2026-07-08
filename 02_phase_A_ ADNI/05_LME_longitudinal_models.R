# Longitudinal Mixed Effects Models

library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(pROC)

select <- dplyr::select

cat("Participants:", length(unique(dat_clean$RID)), "\n")
cat("Observations:", nrow(dat_clean), "\n")

# LONGITUDINAL DATASET

visit_counts <- dat_clean %>%
  group_by(RID) %>%
  summarise(n_visits = n(), .groups = "drop")

multi_visit <- visit_counts %>%
  filter(n_visits >= 2)

dat_long <- dat_clean %>%
  filter(RID %in% multi_visit$RID)

cat("Participants with 2+ visits:", length(unique(dat_long$RID)), "\n")
cat("Total observations:", nrow(dat_long), "\n")

dat_long <- dat_long %>%
  mutate(
    SEX_num   = as.numeric(as.factor(PTGENDER)),
    APOE4_num = as.numeric(as.factor(APOE4))
  )

# PRIMARY MODEL — ENTORHINAL CORTEX

dat_model1 <- dat_long %>%
  select(RID, TIME_YEARS, EC_bilateral,
         CSF_PTAU_bl, AGE, SEX_num,
         PTEDUCAT, APOE4_num) %>%
  na.omit()

cat("Participants in primary model:", length(unique(dat_model1$RID)), "\n")
cat("Observations in primary model:", nrow(dat_model1), "\n")

model_EC_ptau <- lmer(
  EC_bilateral ~
    TIME_YEARS * CSF_PTAU_bl +
    AGE + SEX_num + PTEDUCAT + APOE4_num +
    (1 + TIME_YEARS | RID),
  data    = dat_model1,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa"))

fixed_effects <- as.data.frame(coef(summary(model_EC_ptau)))

cat("=== FIXED EFFECTS TABLE (Entorhinal Cortex) ===\n")
print(round(fixed_effects, 4))

dir.create(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs",
  showWarnings = FALSE,
  recursive    = TRUE)

write.csv(fixed_effects,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table2_LME_EC_ptau.csv")

cat("EC model results saved\n")

# HIPPOCAMPAL MODEL

dat_model2 <- dat_long %>%
  select(RID, TIME_YEARS, HIPP_total,
         CSF_PTAU_bl, AGE, SEX_num,
         PTEDUCAT, APOE4_num) %>%
  na.omit()

cat("Participants in hippocampal model:", length(unique(dat_model2$RID)), "\n")

model_HIPP_ptau <- lmer(
  HIPP_total ~
    TIME_YEARS * CSF_PTAU_bl +
    AGE + SEX_num + PTEDUCAT + APOE4_num +
    (1 + TIME_YEARS | RID),
  data    = dat_model2,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa"))

fixed_hipp <- as.data.frame(coef(summary(model_HIPP_ptau)))

cat("\n=== HIPPOCAMPAL MODEL RESULTS ===\n")
print(round(fixed_hipp, 4))

write.csv(fixed_hipp,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table2_LME_HIPP_ptau.csv")

cat("Hippocampal model saved\n")

# FIGURE 2 — EC TRAJECTORY BY p-tau TERTILE

dat_plot <- dat_long %>%
  filter(!is.na(CSF_PTAU_bl),
         !is.na(EC_bilateral)) %>%
  mutate(
    ptau_tertile = ntile(CSF_PTAU_bl, 3),
    ptau_group   = factor(ptau_tertile,
                          labels = c("Low p-tau181",
                                     "Medium p-tau181",
                                     "High p-tau181")))

cat("Participants per tertile group:\n")
print(table(dat_plot$ptau_group))

traj <- dat_plot %>%
  mutate(TIME_round = round(TIME_YEARS, 0)) %>%
  group_by(ptau_group, TIME_round) %>%
  summarise(
    mean_EC = mean(EC_bilateral, na.rm = TRUE),
    se_EC   = sd(EC_bilateral, na.rm = TRUE) / sqrt(n()),
    n       = n(),
    .groups = "drop") %>%
  filter(n >= 3)

cat("\nTrajectory data points:\n")
print(traj)

fig2 <- ggplot(traj,
               aes(x     = TIME_round,
                   y     = mean_EC,
                   color = ptau_group,
                   fill  = ptau_group)) +
  geom_ribbon(
    aes(ymin = mean_EC - se_EC,
        ymax = mean_EC + se_EC),
    alpha = 0.2,
    color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c("#2196F3", "#FF9800", "#F44336")) +
  scale_fill_manual(
    values = c("#2196F3", "#FF9800", "#F44336")) +
  labs(
    title    = "Entorhinal Cortex Trajectories by Baseline CSF p-tau181",
    subtitle = "Mean \u00B1 SE | ADNI longitudinal sample",
    x        = "Years from Baseline",
    y        = "Entorhinal Cortex Thickness (mm)",
    color    = "Baseline p-tau181",
    fill     = "Baseline p-tau181") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

dir.create(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures",
  showWarnings = FALSE,
  recursive    = TRUE)

ggsave(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig2_EC_trajectory.pdf",
  fig2, width = 8, height = 6, dpi = 300)

ggsave(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig2_EC_trajectory.png",
  fig2, width = 8, height = 6, dpi = 300)

cat("Figure 2 saved\n")

# ROC ANALYSIS

ec_change <- dat_clean %>%
  filter(!is.na(EC_bilateral)) %>%
  arrange(RID, TIME_YEARS) %>%
  group_by(RID) %>%
  summarise(
    EC_baseline  = first(EC_bilateral),
    EC_followup  = EC_bilateral[which.min(abs(TIME_YEARS - 2))],
    time_at_followup = TIME_YEARS[which.min(abs(TIME_YEARS - 2))],
    has_followup = any(abs(TIME_YEARS - 2) < 0.75),
    .groups = "drop") %>%
  filter(has_followup) %>%
  mutate(EC_change = EC_followup - EC_baseline)

cat("Participants with ~2yr follow-up:", nrow(ec_change), "\n")
cat("\nEC change summary (negative = atrophy):\n")
print(summary(ec_change$EC_change))

threshold <- quantile(ec_change$EC_change, 0.25, na.rm = TRUE)

cat("Atrophy threshold (25th percentile):", round(threshold, 4), "mm\n")

ec_change <- ec_change %>%
  mutate(atrophy = as.integer(EC_change <= threshold))

cat("Atrophy cases:", sum(ec_change$atrophy), "\n")
cat("Stable cases:", sum(ec_change$atrophy == 0), "\n")

dat_bl_bio <- dat_clean %>%
  group_by(RID) %>%
  slice(1) %>%
  ungroup() %>%
  select(RID, CSF_PTAU_bl, pT217_bl, AB4240_bl, AGE, DX_bl)

roc_dat_csf <- ec_change %>%
  left_join(dat_bl_bio, by = "RID") %>%
  filter(!is.na(CSF_PTAU_bl))

cat("CSF p-tau181 ROC sample:", nrow(roc_dat_csf), "\n")
cat("Atrophy cases:", sum(roc_dat_csf$atrophy), "\n")

roc_csf <- roc(
  roc_dat_csf$atrophy,
  roc_dat_csf$CSF_PTAU_bl,
  ci        = TRUE,
  ci.method = "bootstrap",
  boot.n    = 2000,
  quiet     = TRUE)

auc_value <- as.numeric(auc(roc_csf))
ci_values <- ci(roc_csf)

cat("=== CSF p-tau181 ROC RESULTS ===\n")
cat("N:", nrow(roc_dat_csf), "\n")
cat("Atrophy cases:", sum(roc_dat_csf$atrophy), "\n")
cat("AUC:", round(auc_value, 3), "\n")
cat("95% CI: [", round(ci_values[1], 3), ",", round(ci_values[3], 3), "]\n")

best_coords <- coords(roc_csf, "best",
                      best.method = "youden",
                      ret = c("threshold", "sensitivity", "specificity"))

cat("\nOptimal threshold:", round(best_coords$threshold, 2), "\n")
cat("Sensitivity:", round(best_coords$sensitivity, 3), "\n")
cat("Specificity:", round(best_coords$specificity, 3), "\n")

table3 <- data.frame(
  Biomarker    = "CSF p-tau181",
  N            = nrow(roc_dat_csf),
  Atrophy_n    = sum(roc_dat_csf$atrophy),
  AUC          = round(auc_value, 3),
  CI_lower     = round(ci_values[1], 3),
  CI_upper     = round(ci_values[3], 3),
  Threshold    = round(best_coords$threshold, 2),
  Sensitivity  = round(best_coords$sensitivity, 3),
  Specificity  = round(best_coords$specificity, 3))

write_csv(table3,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table3_ROC_results.csv")

cat("Table 3 saved\n")

# FIGURE 3 — ROC CURVE

roc_plot_data <- data.frame(
  FPR = 1 - roc_csf$specificities,
  TPR = roc_csf$sensitivities)

fig3 <- ggplot(roc_plot_data,
               aes(x = FPR, y = TPR)) +
  geom_line(color = "#1B3A6B", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  labs(
    title    = "ROC Curve: CSF p-tau181 Predicting Entorhinal Atrophy",
    subtitle = paste0("AUC = ", round(auc_value, 3),
                      " [95% CI: ", round(ci_values[1], 3),
                      "-", round(ci_values[3], 3), "] | N = ",
                      nrow(roc_dat_csf)),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)") +
  theme_minimal(base_size = 13)

ggsave(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig3_ROC_curve.pdf",
  fig3, width = 7, height = 6, dpi = 300)

ggsave(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig3_ROC_curve.png",
  fig3, width = 7, height = 6, dpi = 300)

cat("Figure 3 saved\n")

save.image(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/01_data_cleaning/outputs/phase_A_workspace.RData")

cat("Workspace saved\n")

