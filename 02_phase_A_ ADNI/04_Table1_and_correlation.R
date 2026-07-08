# Table 1 and cross-sectioal correlations

library(dplyr)
library(readr)
library(ggplot2)
library(tableone)
library(ppcor)

select <- dplyr::select

cat("Participants:", length(unique(dat_clean$RID)), "\n")
cat("Observations:", nrow(dat_clean), "\n")

# BASELINE DATASET FOR TABLE 1

dat_baseline <- dat_clean %>%
  filter(VISCODE == "bl")

cat("Participants at baseline:", nrow(dat_baseline), "\n")

missing_bl <- dat_clean %>%
  filter(!RID %in% dat_baseline$RID) %>%
  group_by(RID) %>%
  slice(1) %>%
  ungroup()

cat("Participants without bl visit:", nrow(missing_bl), "\n")

dat_table1 <- bind_rows(dat_baseline, missing_bl)

cat("Total for Table 1:", nrow(dat_table1), "\n")

# TABLE 1

tab1_vars <- c(
  "AGE",
  "PTGENDER",
  "PTEDUCAT",
  "APOE4",
  "EC_bilateral",
  "HIPP_total",
  "CSF_PTAU_bl",
  "pT217_bl",
  "AB4240_bl"
)

tab1_cat <- c("PTGENDER", "APOE4")

table1 <- CreateTableOne(
  vars       = tab1_vars,
  strata     = "DX_bl",
  data       = dat_table1,
  factorVars = tab1_cat,
  test       = TRUE,
  addOverall = TRUE)

print(table1,
      showAllLevels = TRUE,
      quote         = FALSE,
      noSpaces      = TRUE,
      nonnormal     = c("CSF_PTAU_bl", "pT217_bl", "AB4240_bl"))


# SAVE TABLE 1

dir.create(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs",
  showWarnings = FALSE,
  recursive    = TRUE)

tab1_export <- print(table1,
                     showAllLevels = TRUE,
                     printToggle   = FALSE,
                     nonnormal     = c("CSF_PTAU_bl", "pT217_bl", "AB4240_bl"))

write.csv(tab1_export,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/Table1_demographics.csv")

cat("Table 1 saved\n")

# CROSS-SECTIONAL PARTIAL CORRELATIONS

dat_corr <- dat_table1 %>%
  mutate(
    SEX_num   = as.numeric(as.factor(PTGENDER)),
    APOE4_num = as.numeric(as.factor(APOE4))
  )

neuro_vars <- c("EC_bilateral", "HIPP_total")

bio_vars <- c("CSF_PTAU_bl",
              "CSF_TAU_bl",
              "CSF_A4240_bl",
              "pT217_bl",
              "AB4240_bl")

covariates <- c("AGE", "SEX_num", "PTEDUCAT", "APOE4_num")

results <- data.frame()

for(dx_group in c("CN","SMC","EMCI","LMCI","AD")){
  
  dat_sub <- dat_corr %>%
    filter(DX_bl == dx_group)
  
  if(nrow(dat_sub) < 20) next
  
  for(nv in neuro_vars){
    for(bv in bio_vars){
      
      mat <- dat_sub %>%
        select(all_of(c(nv, bv, covariates))) %>%
        na.omit()
      
      if(nrow(mat) < 10) next
      
      out <- tryCatch(
        pcor.test(mat[[1]], mat[[2]], mat[, covariates]),
        error = function(e) NULL
      )
      
      if(is.null(out)) next
      
      results <- rbind(results, data.frame(
        diagnosis = dx_group,
        neuro     = nv,
        biomarker = bv,
        r         = round(out$estimate, 3),
        p         = out$p.value,
        n         = nrow(mat),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("Correlation analyses complete\n")
cat("Number of correlations computed:", nrow(results), "\n")

results <- results %>%
  group_by(diagnosis) %>%
  mutate(
    p_fdr = p.adjust(p, method = "BH"),
    sig   = p_fdr < 0.05,
    stars = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      TRUE          ~ ""
    )
  ) %>%
  ungroup()

cat("\n=== PARTIAL CORRELATION RESULTS ===\n")
print(results %>%
        select(diagnosis, neuro, biomarker, r, p, p_fdr, sig) %>%
        arrange(diagnosis, neuro, p_fdr),
      n = 100)

write_csv(results,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/02_phase_A_ADNI/outputs/TableS1_correlations.csv")

cat("\nResults saved\n")


#  FIGURE 1; CORRELATION HEATMAP

results_plot <- results %>%
  mutate(
    neuro_label = recode(neuro,
                         "EC_bilateral" = "Entorhinal Cx",
                         "HIPP_total"   = "Hippocampus"),
    bio_label = recode(biomarker,
                       "CSF_PTAU_bl"  = "CSF p-tau181",
                       "CSF_TAU_bl"   = "CSF t-tau",
                       "CSF_A4240_bl" = "CSF Abeta42/40",
                       "pT217_bl"     = "Plasma p-tau217",
                       "AB4240_bl"    = "Plasma Abeta42/40"),
    diagnosis = factor(diagnosis,
                       levels = c("CN","SMC","EMCI","LMCI","AD"))
  )

fig1 <- ggplot(results_plot,
               aes(x     = bio_label,
                   y     = neuro_label,
                   fill  = r,
                   label = stars)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low      = "#2196F3",
    mid      = "white",
    high     = "#F44336",
    midpoint = 0,
    limits   = c(-0.6, 0.6),
    name     = "Partial r") +
  facet_wrap(~diagnosis, nrow = 1) +
  labs(
    title    = "Partial Correlations: Brain Regions vs Biomarkers",
    subtitle = "Controlling for age, sex, education, APOE4 | *FDR<0.05",
    x        = "Biomarker",
    y        = "Brain Region") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    strip.background = element_rect(fill = "#1B3A6B"),
    strip.text = element_text(color = "white", face = "bold"))

dir.create(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures",
  showWarnings = FALSE,
  recursive    = TRUE)

ggsave(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig1_correlation_heatmap.pdf",
  fig1, width = 12, height = 5, dpi = 300)

ggsave(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/05_figures/Fig1_correlation_heatmap.png",
  fig1, width = 12, height = 5, dpi = 300)

cat("Figure 1 saved\n")

save.image(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/01_data_cleaning/outputs/phase_A_workspace.RData")

cat("Workspace saved\n")
