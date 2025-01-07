# ==============================================================================
# Script Title: 01_exploratory_analysis.r
# Description: This script performs an exploratory analysis of the pseq object
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
set.seed(123)
options(warn = -1)

# 0.2 Load Packages
# Uncomment and add your packages here
library(MicrobiomeStat)
library(tibble)
library(dplyr)
library(ggpubr)
library(vegan)
library(microbiome)
library(ggfortify)

# 0.3 Declare or Load Custom Functions
# source("")

# 0.4 Declare Variables
# Inputs
input_path <- "02_preprocess/data/"

# Outputs

# Arguments
level <- "Genus"
experiment <- "MS"
levels <- c("MS", "No_MS")

# ==============================================================================
# 1. Load phyloseq object and select columns of interest
# ==============================================================================
# 1.1 Load pseq object
phy <- readRDS(file = paste0(input_path, "phy_",
                             level, ".rds"))

# 1.2 Select columns of interest
df <-  as.data.frame(as.matrix(sample_data(phy), rownames = NA))
df <- df %>%
  select(MS, GS, BMI, DT2_P_H, Aff_H)
sample_data(phy) <- df

# 1.3 Change name of the column experiment as Status to perform the analysis
colnames(phy@sam_data)[colnames(phy@sam_data) == experiment] <- "Status"
phy@sam_data$Status <- factor(x = as.factor(phy@sam_data$Status),
                             levels = levels)
table(phy@sam_data$Status)
target <- "Status"

# ==============================================================================
# 2. Perform Alpha Diversity Analysis
# ==============================================================================
# 2.1 Alpha Diversity Analysis using Shannon and Simpson
rich <- estimate_richness(phy)
pwt <- pairwise.wilcox.test(rich$Shannon, p.adjust.method = "none",
                            sample_data(phy)$Status)
pwt2 <- pairwise.wilcox.test(rich$Simpson, p.adjust.method = "none",
                             sample_data(phy)$Status)
pwt
pwt2

# ==============================================================================
# 3. Perform Beta Diversity Analysis
# ==============================================================================
# 3.1 Beta diversity analysis using Bray-Curtis 
dist <- phyloseq::distance(phy, method = "bray", weighted = FALSE)
adonis2(dist ~ sample_data(phy)$Status)

# 3.2 Beta diversity analysis using Chao
dist2 <- phyloseq::distance(phy, method = "chao", weighted = FALSE)
adonis2(dist2 ~ sample_data(phy)$Status)

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()