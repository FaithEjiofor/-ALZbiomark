setwd("/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/00_data_downloads")

library(dplyr)
library(readr)
library(tidyr)
select <- dplyr::select

cat("Checking all data objects are loaded:\n")
cat("adnimerge loaded:", exists("adnimerge"), "\n")
cat("fs_long loaded:  ", exists("fs_long"),   "\n")
cat("csf loaded:      ", exists("csf"),       "\n")
cat("plasma loaded:   ", exists("plasma"),    "\n")
cat("apoe loaded:     ", exists("apoe"),      "\n")
cat("dx loaded:       ", exists("dx"),        "\n")
cat("mmse loaded:     ", exists("mmse"),      "\n")

# -----------------------------------------------------------
# SECTION 2: SELECT COLUMNS
# -----------------------------------------------------------

adnimerge_select <- adnimerge %>%
  select(RID, COLPROT, PTID, VISCODE, EXAMDATE,
         DX_bl, AGE, PTGENDER, PTEDUCAT, APOE4,
         ABETA, FDG) %>%
  mutate(EXAMDATE = as.Date(EXAMDATE))

cat("ADNIMERGE selected columns:", ncol(adnimerge_select), "\n")
cat("ADNIMERGE rows:", nrow(adnimerge_select), "\n")

fs_select <- fs_long %>%
  select(RID, VISCODE, VISCODE2, EXAMDATE,
         FLDSTRENG,
         ST24TA,
         ST83TA,
         ST29SV,
         ST88SV,
         ST1SV)

cat("FreeSurfer selected columns:", ncol(fs_select), "\n")
cat("FreeSurfer rows:", nrow(fs_select), "\n")

csf_select <- csf %>%
  select(RID, VISCODE2, EXAMDATE,
         ABETA,
         TAU,
         PTAU,
         A4240)

cat("CSF selected columns:", ncol(csf_select), "\n")
cat("CSF rows:", nrow(csf_select), "\n")

plasma_select <- plasma %>%
  select(RID, VISCODE2, EXAMDATE,
         pT217_C2N,
         AB42_C2N,
         AB40_C2N,
         AB42_AB40_C2N)

cat("Plasma selected columns:", ncol(plasma_select), "\n")
cat("Plasma rows:", nrow(plasma_select), "\n")

apoe_select <- apoe %>%
  select(RID, GENOTYPE)

cat("APOE selected columns:", ncol(apoe_select), "\n")
cat("APOE rows:", nrow(apoe_select), "\n")

# -----------------------------------------------------------
# SECTION 3: DERIVED VARIABLES
# -----------------------------------------------------------

fs_select <- fs_select %>%
  mutate(
    EC_bilateral = rowMeans(cbind(ST24TA, ST83TA), na.rm = TRUE),
    HIPP_total   = rowSums(cbind(ST29SV, ST88SV), na.rm = TRUE),
    EXAMDATE     = as.Date(EXAMDATE)
  )

cat("Bilateral EC column created:", "EC_bilateral" %in% names(fs_select), "\n")
cat("Total hippocampus column created:", "HIPP_total" %in% names(fs_select), "\n")

apoe_select <- apoe_select %>%
  mutate(
    APOE4_derived = case_when(
      grepl("4", GENOTYPE) ~ 1,
      !grepl("4", GENOTYPE) ~ 0,
      TRUE ~ NA_real_
    )
  )

cat("\nAPOE4 carrier distribution from GENOTYPE:\n")
print(table(apoe_select$APOE4_derived, useNA = "always"))

apoe_clean <- apoe_select %>%
  select(RID, APOE4_derived) %>%
  distinct(RID, .keep_all = TRUE)

csf_select <- csf_select %>%
  mutate(
    TAU_log   = ifelse(TAU > 0, log(TAU), NA_real_),
    PTAU_log  = ifelse(PTAU > 0, log(PTAU), NA_real_),
    ABETA_log = ifelse(ABETA > 0, log(ABETA), NA_real_)
  )

plasma_select <- plasma_select %>%
  mutate(
    pT217_log = ifelse(pT217_C2N > 0, log(pT217_C2N), NA_real_)
  )

cat("Log transformations created successfully\n")

# -----------------------------------------------------------
# SECTION 4: MERGE ALL FILES
# -----------------------------------------------------------

dat_step1 <- adnimerge_select %>%
  left_join(
    fs_select %>%
      select(RID, VISCODE, FLDSTRENG,
             EC_bilateral, HIPP_total, ST1SV),
    by = c("RID", "VISCODE")
  )

cat("Step 1 complete:\n")
cat("Rows:", nrow(dat_step1), "\n")
cat("EC_bilateral non-missing:",
    sum(!is.na(dat_step1$EC_bilateral)), "\n")

csf_baseline <- csf_select %>%
  arrange(RID, VISCODE2) %>%
  group_by(RID) %>%
  slice(1) %>%
  ungroup() %>%
  select(RID,
         CSF_ABETA_bl    = ABETA,
         CSF_TAU_bl      = TAU,
         CSF_PTAU_bl     = PTAU,
         CSF_A4240_bl    = A4240,
         CSF_TAU_log_bl  = TAU_log,
         CSF_PTAU_log_bl = PTAU_log)

cat("CSF baseline values extracted:\n")
cat("Participants with CSF data:", nrow(csf_baseline), "\n")
cat("CSF PTAU non-missing:",
    sum(!is.na(csf_baseline$CSF_PTAU_bl)), "\n")

plasma_baseline <- plasma_select %>%
  arrange(RID, VISCODE2) %>%
  group_by(RID) %>%
  slice(1) %>%
  ungroup() %>%
  select(RID,
         pT217_bl  = pT217_C2N,
         AB42_bl   = AB42_C2N,
         AB40_bl   = AB40_C2N,
         AB4240_bl = AB42_AB40_C2N)

cat("Plasma baseline values extracted:\n")
cat("Participants with plasma data:", nrow(plasma_baseline), "\n")
cat("Plasma pT217 non-missing:",
    sum(!is.na(plasma_baseline$pT217_bl)), "\n")

dat_step2 <- dat_step1 %>%
  left_join(csf_baseline,    by = "RID") %>%
  left_join(plasma_baseline, by = "RID")

cat("After biomarker merge:\n")
cat("Rows:", nrow(dat_step2), "\n")
cat("Unique participants:", length(unique(dat_step2$RID)), "\n")
cat("CSF PTAU non-missing:", sum(!is.na(dat_step2$CSF_PTAU_bl)), "\n")
cat("Plasma pT217 non-missing:", sum(!is.na(dat_step2$pT217_bl)), "\n")

dat_merged <- dat_step2 %>%
  left_join(
    apoe_clean %>% select(RID, APOE4_derived),
    by = "RID"
  )

cat("Final merged dataset:\n")
cat("Rows:", nrow(dat_merged), "\n")
cat("Unique participants:", length(unique(dat_merged$RID)), "\n")
cat("Columns:", ncol(dat_merged), "\n")

# -----------------------------------------------------------
# SECTION 5 & 6: TIME VARIABLE + INCLUSION CRITERIA
# (built directly from dat_merged in one pipeline)
# -----------------------------------------------------------

dat_clean <- dat_merged %>%
  mutate(EXAMDATE = as.Date(EXAMDATE)) %>%
  group_by(RID) %>%
  mutate(
    BASELINE_DATE = min(EXAMDATE, na.rm = TRUE),
    TIME_YEARS = as.numeric(EXAMDATE - BASELINE_DATE) / 365.25
  ) %>%
  ungroup() %>%
  filter(AGE >= 40) %>%
  filter(DX_bl %in% c("CN", "SMC", "EMCI", "LMCI", "AD")) %>%
  filter(!is.na(EC_bilateral))

cat("============================================\n")
cat("FINAL CLEAN DATASET\n")
cat("============================================\n")
cat("Participants:", length(unique(dat_clean$RID)), "\n")
cat("Observations:", nrow(dat_clean), "\n")
cat("Columns:", ncol(dat_clean), "\n")
cat("============================================\n")

# -----------------------------------------------------------
# SECTION 7: SAVE
# -----------------------------------------------------------

dir.create(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/01_data_cleaning/outputs",
  showWarnings = FALSE,
  recursive    = TRUE)

write_csv(dat_clean,
          "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/01_data_cleaning/outputs/adni_clean_analysis.csv")

save.image(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/01_data_cleaning/outputs/phase_A_workspace.RData")

cat("Clean dataset and workspace saved successfully\n")
save.image(
  "/home/faithejiofor15/alzheimers-biomarker-neuroanatomy/01_data_cleaning/outputs/phase_A_workspace.RData")

cat("Clean dataset and workspace saved successfully\n")


