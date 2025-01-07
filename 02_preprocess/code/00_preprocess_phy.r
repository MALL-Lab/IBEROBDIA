# ==============================================================================
# Script Title: 00_preprocess_phy.r
# Description: This script preprocesses the phyloseq object
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
library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)

# 0.3 Declare or Load Custom Functions
# source("")

# 0.4 Declare Variables
# Arguments
tl <- "TL_251" # Experiment name
rank <- 7 # Taxonomic rank to agglomerate

# Inputs
input_path <- paste0("01_sequencing_data/data/", TL, "/", TL, "_phy.rds")

# Outputs
out_path <- ("02_preprocess/data/phy")
# ==============================================================================
# 1. Load phyloseq object and summarize
# ==============================================================================
# 1.1 Load pseq object
phy <- readRDS(input_path)

# 1.2 Summary
ntaxa(phy)
nsamples(phy)
sample_names(phy)[1:5]
rank_names(phy)
sample_variables(phy)
otu_table(phy)[1:5, 1:5]
tax_table(phy)[1:5, 1:7]

# ==============================================================================
# 2. Agglomerate by taxonomic rank
# ==============================================================================
g_phy <- tax_glom(physeq = phy, taxrank = rank_names(phy)[rank])

# ==============================================================================
# 3. Filter out ASVs that do not meet certain criteria.
# ==============================================================================
# Remove taxa that are not seen more than 3 times in at least 20% of the samples
gp <- filter_taxa(g_phy, function(x) sum(x > 3) > (0.2 * length(x)), TRUE)

# ==============================================================================
# 4.  Remove samples that do not have a minimum of 20 counts
# ==============================================================================
gp <- prune_samples(sample_sums(gp) >= 20, gp)

# ==============================================================================
# 5. Check and save the phyloseq object
# ==============================================================================
any(taxa_sums(phy) == 0)
saveRDS(object = gp,
        file = paste0(input_path, "_", rank_names(phy)[rank], ".rds"))

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()