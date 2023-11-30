# Script to make DEA analysis 
# -------------------------------------
# 3 variables:
# level :"Genus" or "Species"
#
# experiment:"DT2_All" "DT2_NoDT2 "DT2_P_H" "DT2_H" "Aff_H" "MS" "GS" "BMI" "IR"
# "DPH_MS" "DPH_GS" "DPH_BMI" "DPH_GS_BMI" "AffH_MS" "AffH_GS" "AffH_BMI"
# "AffH_GS_BMI" "AffH_MS_GS_BMI" "MS_GS" "MS_BMI" "MS_IR" "MS_GS_BMI" "BMI_IR"

require(MicrobiomeStat)
require(tibble)
require(dplyr)
library(ggpubr)
library(vegan)
library(microbiome)
library(ggfortify)
set.seed(123)

# Declare variables 
level <- "Genus"
experiment <- "DT2_P_H"
levels <- c("Healthy", "PreDT2","DT2" )
# Load phyloseq
phy <- readRDS(file = paste0("02_preprocess/data/phy_",
                             level,"_v2.rds"))
df <-  as.data.frame(as.matrix(sample_data(phy), rownames = NA))
df <- df %>% 
  select(MS, GS, BMI, DT2_P_H, Aff_H)
sample_data(phy) <- df

# Change name of the column experiment as Status to perform the analysis
colnames(phy@sam_data)[colnames(phy@sam_data) == experiment] <- "Status"
phy@sam_data$Status <-factor(x = as.factor(phy@sam_data$Status),
                             levels = levels)
table(phy@sam_data$Status)
target <- "Status"
 
# Alpha Diversity 
rich <- estimate_richness(phy)
pwt <- pairwise.wilcox.test(rich$Shannon, p.adjust.method = "none",
                            sample_data(phy)$Status)
pwt2 <- pairwise.wilcox.test(rich$Simpson, p.adjust.method = "none",
                            sample_data(phy)$Status)
pwt
pwt2

# Beta diversity
dist = phyloseq::distance(phy,
                          method ="bray",
                          weighted = F)
dist2 = phyloseq::distance(phy,
                          method ="chao",
                          weighted = F)
adonis2(dist ~ sample_data(phy)$Status)
adonis2(dist2 ~ sample_data(phy)$Status)
