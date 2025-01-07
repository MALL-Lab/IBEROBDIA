# ==============================================================================
# Script Title: 00_preprocess_metadata
# Description: This script preprocesses the metadata of the samples.
# ==============================================================================
# Author/s: @DiegoFE94 (GitHub)
# Affiliation: Machine Learning for Life Sciences Laboratory (MALL)
# Email: diego.fedreira@udc.es
# Date Created: 2024 (Last Update on January 2025)
# ==============================================================================
# 0. Set Up
# ==============================================================================
# 0.1 Clean Environment
rm(list = ls())
set.seed(1965)

# 0.2 Load Packages
require(tibble)
require(readxl)
require(dplyr)
library(stringr)
library(data.table)

# 0.3 Declare or Load Custom Functions
# source("")

# 0.4 Declare Variables

# Inputs
path_samples <- "extdata/FASTQ/"
path_clinical_data <- "extdata/clinical_data.xlsx"

# Outputs
path_metadata <- "00_preprocess_metadata/data/clinical_data.rds"

# Arguments

# ==============================================================================
# 1. Clinical Data Extraction
# ==============================================================================
# 1.1 Extract sample names
names <- list.files(path = path_samples, pattern = ".fastq.gz")
vars <- list()
for (i in seq_along(names)) {
  vars[i] <- str_remove(names[i], ".fastq.gz")
  vars[i] <- str_split(vars[i], "_", n = Inf, simplify = FALSE)
}

# 1.2 Extract clinical data
metalist <- list()
for (i in seq_along(vars)) {
  if (length(vars[[i]]) == 3) {
    metalist[[i]] <- c(vars[[i]][1], vars[[i]][2], "Healthy", vars[[i]][3])
  }else {
    metalist[[i]] <- vars[[i]]
  }
}

# 1.3 Create clinical dataframe
clinical <- as.data.frame(do.call(rbind, metalist))
colnames(clinical) <- c("Sample_ID", "CI", "DT2_All", "Gender")
clinical <-  as_tibble(x = clinical, rownames = NA)

# 1.4 Preprocess clinical data
clinical <- clinical %>%
  select(Sample_ID, CI, Gender, DT2_All) %>%
  mutate(DT2_P_H = ifelse(DT2_All == "GAA", "PreDT2",
                          ifelse(DT2_All == "TGA", "PreDT2",
                                 DT2_All))) %>%
  mutate(Aff_H = ifelse(DT2_All == "GAA", "Affected",
                        ifelse(DT2_All == "TGA", "Affected",
                               ifelse(DT2_All == "DT2", "Affected",
                                      DT2_All))))

# ==============================================================================
# 2. Merge with the rest of the clinical data
# ==============================================================================
# 2.1 Load the rest of the clinical data
clinical_v2 <- readxl::read_xlsx(path_clinical_data)
clinical_v2 <- as_tibble(clinical_v2)
clinical_v2 <- clinical_v2 %>%
  select(Sample_ID, MS, GS, BMI, IR)

# 2.2 Merge the clinical data
full_df <- dplyr::full_join(clinical, clinical_v2, by = "Sample_ID")
# ==============================================================================
# 3. Process merge data
# ==============================================================================
# 3.1 Preprocess the metadata
metadata <- full_df %>%
  # Delete 084CMLES sample (No FASTQ)
  filter(Sample_ID != "084CMLES") %>%
  # Change Metabolic Syndrome variable
  mutate(MS = ifelse(MS == "No", "No_MS", MS)) %>%
  mutate(MS = ifelse(MS == "Sí", "MS", MS)) %>%
  # Change Status of Glucosevariable
  mutate(GS = ifelse(GS == "No", "NG", GS)) %>%
  mutate(GS = ifelse(GS == "Sí", "AG", GS)) %>%
  # Change Body Mass Index variable
  mutate(BMI = ifelse(BMI == "Sobrepeso 2", "OB", BMI)) %>%
  mutate(BMI = ifelse(BMI == "Obesidad", "OB", BMI)) %>%
  mutate(BMI = ifelse(BMI == "Normopeso", "NW", BMI)) %>%
  # Change Insuline Resistanece variable
  mutate(IR = ifelse(IR == "Resistencia significativa", "SR", IR)) %>%
  mutate(IR = ifelse(IR == "Resistencia temprana", "ER", IR)) %>%
  mutate(IR = ifelse(IR == "Resistencia Temprana", "ER", IR)) %>%
  mutate(IR = ifelse(IR == "Rango Saludable", "No_IR", IR)) %>%
  mutate(IR = ifelse(IR == "Rango saludable", "No_IR", IR)) %>%
  column_to_rownames(var = "Sample_ID")

# ==============================================================================
# 4. Save the clinical data
# ==============================================================================
saveRDS(object = metadata, file = path_metadata)

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()
# Preprocess_metadata
