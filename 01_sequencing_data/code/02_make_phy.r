# ==============================================================================
# Script Title: 02_make_phy.r
# Description: This script creates a phyloseq object
# ==============================================================================
# Author/s: @DiegoFE94 (GitHub)
# Affiliation: Machine Learning for Life Sciences Laboratory (MALL)
# Email: diego.fedreira@udc.es
# Date Created: # Date Created: 2024 (Last Update on January 2025)
# ==============================================================================
# 0. Set Up
# ==============================================================================
# 0.1 Clean Environment
rm(list = ls())
set.seed(1965)

# 0.2 Load Packages
library(phyloseq)
library(dplyr)

# 0.3 Declare or Load Custom Functions
# source("")

# 0.4 Declare Variables
# Inputs
metadata_input <- "00_preprocess_metadata/data/clinical_data.rds"
tables_input <- "01_sequencing_data/data/"

# Outputs

# Arguments

# ==============================================================================
# 1. Load clinical data and prepare paths
# ==============================================================================
# 1.1 Load metadata
metadata <- readRDS(metadata_input)

# 1.2 Prepare paths
truncL <- sort(list.files(tables_input, pattern = "TL_251"))
nams <- paste(truncL, "phy", sep = "_")

# ==============================================================================
# 2. Create phyloseq object
# ==============================================================================
# 2.1 create a list for pseq objects
phy_list <- list()

# 2.2 Pseq objects construction
for (i in seq_along(truncL)) {
  # 2.2.1 Load and merge all data in pseq object
  tax <- readRDS(file = paste(path, truncL[i], "tax_table.rds", sep = "/"))
  otu <- readRDS(file = paste(path, truncL[i], "otu_table.rds", sep = "/"))
  otu <- phyloseq::t(otu)
  ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
                 sample_data(metadata),
                 tax_table(tax))

  # 2.2.2 Shorten the name of our ASVs
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

  # 2.2.3 Correctly rename the species for a correct merge later on.
  ## This is done in order to have a complete species name and
  ## not to make mistakes when agglomerating.
  ## Replace NAs for 0 to filter spss in order to add genus name before
  tax <- data.frame(tax_table(ps))
  tax$Species[is.na(tax$Species)] <- 0
  ## Add genus name to Species names
  tax <- tax %>%
    mutate(Species = ifelse(Species != 0, paste0(Genus, "_", Species), Species))

  ## Replace 0 for NA again
  tax <- tax %>%
    mutate(Species = ifelse(Species == 0, NA, Species))
  tax <- as.matrix(tax)
  tax_table(ps) <- tax
  phy_list[[i]] <- ps
}

# ==============================================================================
# 3. Rename and save the pseq objects
# ==============================================================================
names(phy_list) <- nams
for (i in seq_along(phy_list)) {
  saveRDS(object = phy_list[[i]],
          file = paste(path, truncL[i], paste0(names(phy_list[i]), ".rds"),
                       sep = "/"))
}

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()