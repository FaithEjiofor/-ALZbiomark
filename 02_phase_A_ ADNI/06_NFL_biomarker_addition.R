#PHASE A RERUN: Adding CSF NfL and Plasma NfL_Q (update) 

library(dplyr)
library(readr)
library(ggplot2)
library(tableone)
library(ppcor)
library(lme4)
library(lmerTest)
library(pROC)

select <- dplyr::select


# SECTION 1: LOAD AND INSPECT NEW BIOMARKER FILES


setwd("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/00_data_downloads")

# CSF NfL (Blennow Lab, single platform, no missing-value codes)
csf_nfl_file <- read_csv("BLENNOWCSFNFL_28Jul2026.csv", show_col_types = FALSE)
cat("CSF NfL file loaded:", nrow(csf_nfl_file), "rows,",
    ncol(csf_nfl_file), "columns\n")

# GFAP/NfL plasma file (UPENN Fujirebio & Quanterix platforms)
gfap_file <- read_csv("UPENN_PLASMA_FUJIREBIO_QUANTERIX_28Jul2026.csv",
                      show_col_types = FALSE)
cat("UPENN plasma file loaded:", nrow(gfap_file), "rows,",
    ncol(gfap_file), "columns\n")


# SECTION 2: RECODE MISSING-VALUE CODES (-4, -5) TO NA

# ADNI/UPENN files use -4 and -5 as standard missing/not-
# applicable codes, not real biomarker concentrations. GFAP_F
# and GFAP_Q both required this; NfL_F and NfL_Q on the same
# platform files require identical treatment. CSF NfL (Blennow
# file) was checked and contains no such codes - used as-is.

gfap_file$NfL_F_clean <- ifelse(gfap_file$NfL_F %in% c(-4, -5),
                                NA, gfap_file$NfL_F)
gfap_file$NfL_Q_clean <- ifelse(gfap_file$NfL_Q %in% c(-4, -5),
                                NA, gfap_file$NfL_Q)

cat("NfL_F_clean summary:\n")
print(summary(gfap_file$NfL_F_clean))
cat("\nNfL_Q_clean summary:\n")
print(summary(gfap_file$NfL_Q_clean))


# SECTION 3: VERIFY PLATFORM AGREEMENT (F vs Q)

# _F = Fujirebio Lumipulse G1200 platform
# _Q = Quanterix Simoa HD-X platform
# Per ADNI/UPENN documentation, these are two independent
# assay platforms, not "final vs QC" values. Correlation
# checked as a sanity check before deciding which to use.

nfl_platform_cor <- cor(gfap_file$NfL_F_clean, gfap_file$NfL_Q_clean,
                        use = "complete.obs")
cat("Correlation between NfL_F and NfL_Q:", round(nfl_platform_cor, 3), "\n")
# Result: r = 0.925 - strong agreement between platforms.
# Quanterix (Q) selected as primary plasma NfL variable due to
# substantially larger sample coverage (1727 vs 1247 valid
# values) despite both being legitimate measurements.


# SECTION 4: EXTRACT BASELINE VALUES AND MERGE INTO dat_clean


csf_nfl_baseline <- csf_nfl_file %>%
  arrange(RID, EXAMDATE) %>%
  distinct(RID, .keep_all = TRUE) %>%
  select(RID, CSF_NfL_bl = CSFNFL)

cat("CSF NfL baseline participants:", nrow(csf_nfl_baseline), "\n")

plasma_nfl_q_baseline <- gfap_file %>%
  filter(!is.na(NfL_Q_clean)) %>%
  arrange(RID, EXAMDATE) %>%
  distinct(RID, .keep_all = TRUE) %>%
  select(RID, Plasma_NfL_Q_bl = NfL_Q_clean)

cat("Plasma NfL_Q baseline participants:", nrow(plasma_nfl_q_baseline), "\n")

dat_clean <- dat_clean %>%
  left_join(csf_nfl_baseline, by = "RID") %>%
  left_join(plasma_nfl_q_baseline, by = "RID")

cat("CSF NfL non-missing in dat_clean (unique RID):",
    length(unique(dat_clean$RID[!is.na(dat_clean$CSF_NfL_bl)])), "\n")
cat("Plasma NfL_Q non-missing in dat_clean (unique RID):",
    length(unique(dat_clean$RID[!is.na(dat_clean$Plasma_NfL_Q_bl)])), "\n")
# Confirmed: 383 and 182 respectively, out of 760 total participants.


# SECTION 5: UPDATED TABLE 1


dat_baseline <- dat_clean %>% filter(VISCODE == "bl")

missing_bl <- dat_clean %>%
  filter(!RID %in% dat_baseline$RID) %>%
  arrange(RID, TIME_YEARS) %>%
  group_by(RID) %>%
  slice(1) %>%
  ungroup()

dat_table1 <- bind_rows(dat_baseline, missing_bl)

tab1_vars <- c("AGE", "PTGENDER", "PTEDUCAT", "APOE4",
               "EC_bilateral", "HIPP_total",
               "CSF_PTAU_bl", "pT217_bl", "AB4240_bl",
               "CSF_NfL_bl", "Plasma_NfL_Q_bl")
tab1_cat <- c("PTGENDER", "APOE4")

table1_updated <- CreateTableOne(
  vars = tab1_vars, strata = "DX_bl", data = dat_table1,
  factorVars = tab1_cat, test = TRUE, addOverall = TRUE)

dir.create(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs",
  showWarnings = FALSE, recursive = TRUE)

tab1_export_updated <- print(table1_updated, showAllLevels = TRUE,
                             printToggle = FALSE,
                             nonnormal = c("CSF_PTAU_bl", "pT217_bl", "AB4240_bl",
                                           "CSF_NfL_bl", "Plasma_NfL_Q_bl"))

write.csv(tab1_export_updated,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table1_demographics_v2.csv")


# SECTION 6: UPDATED CROSS-SECTIONAL CORRELATIONS


dat_corr <- dat_table1 %>%
  mutate(SEX_num = as.numeric(as.factor(PTGENDER)),
         APOE4_num = as.numeric(as.factor(APOE4)))

neuro_vars <- c("EC_bilateral", "HIPP_total")
bio_vars <- c("CSF_PTAU_bl", "CSF_TAU_bl", "CSF_A4240_bl",
              "pT217_bl", "AB4240_bl",
              "CSF_NfL_bl", "Plasma_NfL_Q_bl")
covariates <- c("AGE", "SEX_num", "PTEDUCAT", "APOE4_num")

results_v2 <- data.frame()

for(dx_group in c("CN","SMC","EMCI","LMCI","AD")){
  dat_sub <- dat_corr %>% filter(DX_bl == dx_group)
  if(nrow(dat_sub) < 20) next
  for(nv in neuro_vars){
    for(bv in bio_vars){
      mat <- dat_sub %>% select(all_of(c(nv, bv, covariates))) %>% na.omit()
      if(nrow(mat) < 10) next
      out <- tryCatch(pcor.test(mat[[1]], mat[[2]], mat[, covariates]),
                      error = function(e) NULL)
      if(is.null(out)) next
      results_v2 <- rbind(results_v2, data.frame(
        diagnosis = dx_group, neuro = nv, biomarker = bv,
        r = round(out$estimate, 3), p = out$p.value,
        n = nrow(mat), stringsAsFactors = FALSE))
    }
  }
}

results_v2 <- results_v2 %>%
  group_by(diagnosis) %>%
  mutate(p_fdr = p.adjust(p, method = "BH"),
         sig = p_fdr < 0.05,
         stars = case_when(p_fdr < 0.001 ~ "***", p_fdr < 0.01 ~ "**",
                           p_fdr < 0.05 ~ "*", TRUE ~ "")) %>%
  ungroup()

write_csv(results_v2,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/TableS1_correlations_v2.csv")

# KEY FINDING: LMCI group, CSF NfL vs EC_bilateral: r=-0.205, p_fdr=0.0287 (sig)
#              LMCI group, CSF NfL vs HIPP_total:  r=-0.299, p_fdr=0.000427 (sig)
# These are the only two correlations surviving FDR correction
# across the entire Phase A analysis (original 5 + new 2 biomarkers).


# SECTION 7: UPDATED FIGURE 1 - CORRELATION HEATMAP


results_plot_v2 <- results_v2 %>%
  mutate(
    neuro_label = recode(neuro, "EC_bilateral" = "Entorhinal Cx",
                         "HIPP_total" = "Hippocampus"),
    bio_label = recode(biomarker,
                       "CSF_PTAU_bl" = "CSF p-tau181", "CSF_TAU_bl" = "CSF t-tau",
                       "CSF_A4240_bl" = "CSF A\u03B242/40", "pT217_bl" = "Plasma p-tau217",
                       "AB4240_bl" = "Plasma A\u03B242/40", "CSF_NfL_bl" = "CSF NfL",
                       "Plasma_NfL_Q_bl" = "Plasma NfL"),
    diagnosis = factor(diagnosis, levels = c("CN","SMC","EMCI","LMCI","AD")))

fig1_v2 <- ggplot(results_plot_v2,
                  aes(x = bio_label, y = neuro_label, fill = r, label = stars)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#2196F3", mid = "white", high = "#F44336",
                       midpoint = 0, limits = c(-0.6, 0.6), name = "Partial r") +
  facet_wrap(~diagnosis, nrow = 1) +
  labs(title = "Partial Correlations: Brain Regions vs Biomarkers (Updated)",
       subtitle = "Controlling for age, sex, education, APOE4 | *FDR<0.05 | Now includes CSF NfL and Plasma NfL",
       x = "Biomarker", y = "Brain Region") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        strip.background = element_rect(fill = "#1B3A6B"),
        strip.text = element_text(color = "white", face = "bold"))

dir.create("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures",
           showWarnings = FALSE, recursive = TRUE)

ggsave("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig1_correlation_heatmap_v2.pdf",
       fig1_v2, width = 16, height = 5, dpi = 300)
ggsave("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig1_correlation_heatmap_v2.png",
       fig1_v2, width = 16, height = 5, dpi = 300)


# SECTION 8: LONGITUDINAL DATASET FOR LME MODELS


visit_counts <- dat_clean %>% group_by(RID) %>%
  summarise(n_visits = n(), .groups = "drop")
multi_visit <- visit_counts %>% filter(n_visits >= 2)

dat_long <- dat_clean %>%
  filter(RID %in% multi_visit$RID) %>%
  mutate(SEX_num = as.numeric(as.factor(PTGENDER)),
         APOE4_num = as.numeric(as.factor(APOE4)))

cat("Participants with 2+ visits:", length(unique(dat_long$RID)), "\n")
cat("CSF NfL available:", length(unique(dat_long$RID[!is.na(dat_long$CSF_NfL_bl)])), "\n")
cat("Plasma NfL available:", length(unique(dat_long$RID[!is.na(dat_long$Plasma_NfL_Q_bl)])), "\n")
# Confirmed: 696 total, 353 with CSF NfL, 174 with Plasma NfL_Q


# SECTION 9: LME MODEL - ENTORHINAL CORTEX ~ CSF NfL

dat_model_csfnfl <- dat_long %>%
  select(RID, TIME_YEARS, EC_bilateral, CSF_NfL_bl,
         AGE, SEX_num, PTEDUCAT, APOE4_num) %>%
  na.omit()

model_EC_csfnfl <- lmer(
  EC_bilateral ~ TIME_YEARS * CSF_NfL_bl + AGE + SEX_num + PTEDUCAT + APOE4_num +
    (1 + TIME_YEARS | RID),
  data = dat_model_csfnfl, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa"))

fixed_EC_csfnfl <- as.data.frame(coef(summary(model_EC_csfnfl)))
write.csv(fixed_EC_csfnfl,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table2_LME_EC_CSFNfL.csv")

# RESULT: Main effect CSF_NfL_bl p<0.0001 (***)
#         Interaction TIME_YEARS:CSF_NfL_bl p=0.0192 (*) - SIGNIFICANT
#         Higher baseline CSF NfL predicts faster entorhinal thinning.


# SECTION 10: LME MODEL - HIPPOCAMPUS ~ CSF NfL

dat_model_csfnfl_hipp <- dat_long %>%
  select(RID, TIME_YEARS, HIPP_total, CSF_NfL_bl,
         AGE, SEX_num, PTEDUCAT, APOE4_num) %>%
  na.omit()

model_HIPP_csfnfl <- lmer(
  HIPP_total ~ TIME_YEARS * CSF_NfL_bl + AGE + SEX_num + PTEDUCAT + APOE4_num +
    (1 + TIME_YEARS | RID),
  data = dat_model_csfnfl_hipp, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa"))

fixed_HIPP_csfnfl <- as.data.frame(coef(summary(model_HIPP_csfnfl)))
write.csv(fixed_HIPP_csfnfl,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table2_LME_HIPP_CSFNfL.csv")

# RESULT: Main effect CSF_NfL_bl p<0.0001 (***)
#         Interaction TIME_YEARS:CSF_NfL_bl p=0.0034 (**) - SIGNIFICANT
#         Strongest longitudinal finding in the entire Phase A analysis.


# SECTION 11: LME MODEL - ENTORHINAL CORTEX ~ PLASMA NfL

dat_model_plasmanfl <- dat_long %>%
  select(RID, TIME_YEARS, EC_bilateral, Plasma_NfL_Q_bl,
         AGE, SEX_num, PTEDUCAT, APOE4_num) %>%
  na.omit()

model_EC_plasmanfl <- lmer(
  EC_bilateral ~ TIME_YEARS * Plasma_NfL_Q_bl + AGE + SEX_num + PTEDUCAT + APOE4_num +
    (1 + TIME_YEARS | RID),
  data = dat_model_plasmanfl, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa"))

fixed_EC_plasmanfl <- as.data.frame(coef(summary(model_EC_plasmanfl)))
write.csv(fixed_EC_plasmanfl,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table2_LME_EC_PlasmaNfL.csv")

# RESULT: Main effect Plasma_NfL_Q_bl p=0.0002 (***)
#         Interaction TIME_YEARS:Plasma_NfL_Q_bl p=0.2919 - NOT significant



# SECTION 12: LME MODEL - HIPPOCAMPUS ~ PLASMA NfL

dat_model_plasmanfl_hipp <- dat_long %>%
  select(RID, TIME_YEARS, HIPP_total, Plasma_NfL_Q_bl,
         AGE, SEX_num, PTEDUCAT, APOE4_num) %>%
  na.omit()

model_HIPP_plasmanfl <- lmer(
  HIPP_total ~ TIME_YEARS * Plasma_NfL_Q_bl + AGE + SEX_num + PTEDUCAT + APOE4_num +
    (1 + TIME_YEARS | RID),
  data = dat_model_plasmanfl_hipp, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa"))

fixed_HIPP_plasmanfl <- as.data.frame(coef(summary(model_HIPP_plasmanfl)))
write.csv(fixed_HIPP_plasmanfl,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table2_LME_HIPP_PlasmaNfL.csv")

# RESULT: Main effect Plasma_NfL_Q_bl p=0.0043 (**)
#         Interaction TIME_YEARS:Plasma_NfL_Q_bl p=0.5002 - NOT significant



# SECTION 13: FIGURE 2 (v2) - EC TRAJECTORY BY CSF NfL TERTILE

dat_plot_nfl <- dat_long %>%
  filter(!is.na(CSF_NfL_bl), !is.na(EC_bilateral)) %>%
  mutate(nfl_tertile = ntile(CSF_NfL_bl, 3),
         nfl_group = factor(nfl_tertile,
                            labels = c("Low CSF NfL", "Medium CSF NfL", "High CSF NfL")))

traj_nfl <- dat_plot_nfl %>%
  mutate(TIME_round = round(TIME_YEARS, 0)) %>%
  group_by(nfl_group, TIME_round) %>%
  summarise(mean_EC = mean(EC_bilateral, na.rm = TRUE),
            se_EC = sd(EC_bilateral, na.rm = TRUE) / sqrt(n()),
            n = n(), .groups = "drop") %>%
  filter(n >= 3)

fig2_nfl <- ggplot(traj_nfl,
                   aes(x = TIME_round, y = mean_EC, color = nfl_group, fill = nfl_group)) +
  geom_ribbon(aes(ymin = mean_EC - se_EC, ymax = mean_EC + se_EC),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
  scale_color_manual(values = c("#2196F3", "#FF9800", "#F44336")) +
  scale_fill_manual(values = c("#2196F3", "#FF9800", "#F44336")) +
  labs(title = "Entorhinal Cortex Trajectories by Baseline CSF NfL",
       subtitle = "Mean \u00B1 SE | ADNI longitudinal sample",
       x = "Years from Baseline", y = "Entorhinal Cortex Thickness (mm)",
       color = "Baseline CSF NfL", fill = "Baseline CSF NfL") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")

ggsave("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig2_EC_trajectory_CSFNfL.pdf",
       fig2_nfl, width = 8, height = 6, dpi = 300)
ggsave("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig2_EC_trajectory_CSFNfL.png",
       fig2_nfl, width = 8, height = 6, dpi = 300)


# SECTION 14: ROC ANALYSIS - CSF NfL AND PLASMA NfL


dat_bl_bio_v2 <- dat_clean %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  select(RID, CSF_NfL_bl, Plasma_NfL_Q_bl, AGE, DX_bl)

roc_dat_csfnfl <- ec_change %>%
  left_join(dat_bl_bio_v2, by = "RID") %>%
  filter(!is.na(CSF_NfL_bl))

roc_csfnfl <- roc(roc_dat_csfnfl$atrophy, roc_dat_csfnfl$CSF_NfL_bl,
                  ci = TRUE, ci.method = "bootstrap", boot.n = 2000, quiet = TRUE)

roc_dat_plasmanfl <- ec_change %>%
  left_join(dat_bl_bio_v2, by = "RID") %>%
  filter(!is.na(Plasma_NfL_Q_bl))

roc_plasmanfl <- roc(roc_dat_plasmanfl$atrophy, roc_dat_plasmanfl$Plasma_NfL_Q_bl,
                     ci = TRUE, ci.method = "bootstrap", boot.n = 2000, quiet = TRUE)

# RESULTS:
# CSF NfL:    N=306, AUC=0.584, 95% CI [0.512, 0.656]
# Plasma NfL: N=155, AUC=0.605, 95% CI [0.510, 0.699]
# (Compare to original CSF p-tau181: N=40, AUC=0.681, 95% CI [0.502, 0.860])
# NOTE: All three lower CI bounds sit near 0.5 - none confidently
# exceed chance-level discrimination for this specific classification task,
# despite CSF NfL's strong longitudinal (LME) findings above. See
# Transparent Changes addendum for discussion of this divergence.

table3_v2 <- data.frame(
  Biomarker = c("CSF p-tau181", "CSF NfL", "Plasma NfL"),
  N = c(40, nrow(roc_dat_csfnfl), nrow(roc_dat_plasmanfl)),
  Atrophy_n = c(13, sum(roc_dat_csfnfl$atrophy), sum(roc_dat_plasmanfl$atrophy)),
  AUC = c(0.681, round(as.numeric(auc(roc_csfnfl)), 3),
          round(as.numeric(auc(roc_plasmanfl)), 3)),
  CI_lower = c(0.502, round(ci(roc_csfnfl)[1], 3), round(ci(roc_plasmanfl)[1], 3)),
  CI_upper = c(0.860, round(ci(roc_csfnfl)[3], 3), round(ci(roc_plasmanfl)[3], 3)))

write_csv(table3_v2,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table3_ROC_results_v2.csv")


# SECTION 15: FIGURE 3 (v2) - ROC CURVES, CSF NfL AND PLASMA NfL

roc_plot_csfnfl <- data.frame(FPR = 1 - roc_csfnfl$specificities,
                              TPR = roc_csfnfl$sensitivities, Biomarker = "CSF NfL")
roc_plot_plasmanfl <- data.frame(FPR = 1 - roc_plasmanfl$specificities,
                                 TPR = roc_plasmanfl$sensitivities, Biomarker = "Plasma NfL")
roc_combined <- bind_rows(roc_plot_csfnfl, roc_plot_plasmanfl)

fig3_nfl <- ggplot(roc_combined, aes(x = FPR, y = TPR, color = Biomarker)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#1B3A6B", "#C0392B")) +
  labs(title = "ROC Curves: CSF NfL and Plasma NfL Predicting Entorhinal Atrophy",
       subtitle = "CSF NfL: AUC=0.584 (N=306) | Plasma NfL: AUC=0.605 (N=155)",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")

ggsave("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig3_ROC_curve_NfL.pdf",
       fig3_nfl, width = 7, height = 6, dpi = 300)
ggsave("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig3_ROC_curve_NfL.png",
       fig3_nfl, width = 7, height = 6, dpi = 300)


# SECTION 16: UPDATED PHASE C ML DATASET 

ml_baseline_v2 <- dat_clean %>%
  filter(VISCODE == "bl" | (!RID %in% dat_clean$RID[dat_clean$VISCODE == "bl"])) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  mutate(SEX_num = as.numeric(as.factor(PTGENDER)),
         APOE4_num = as.numeric(as.factor(APOE4)))

ml_data_v2 <- ml_baseline_v2 %>%
  filter(DX_bl == "CN") %>%
  left_join(conversion_outcome_full, by = "RID") %>%
  select(RID, DX_bl, EC_bilateral, HIPP_total,
         CSF_PTAU_bl, CSF_TAU_bl, CSF_A4240_bl, pT217_bl, AB4240_bl,
         CSF_NfL_bl, Plasma_NfL_Q_bl,
         AGE, SEX_num, PTEDUCAT, APOE4_num, converted, max_followup)

dir.create("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/04_phase_C_ML",
           showWarnings = FALSE, recursive = TRUE)

write_csv(ml_data_v2,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/04_phase_C_ML/ml_baseline_data_v2.csv")

# CONFIRMED: 217 CN participants, 22 conversion events (unchanged from
# original Phase C outcome logic - only predictors were added).
# Model B viability check:
#   CSF NfL:    N=109, 11 events - meets minimum threshold (marginal)
#   Plasma NfL: N=86,  21 events - meets minimum threshold (stronger)


# SECTION 17: SAVE UPDATED WORKSPACE

save.image(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/01_data_cleaning/outputs/phase_A_workspace_v2.RData")

cat("SCRIPT 6 COMPLETE\n")
