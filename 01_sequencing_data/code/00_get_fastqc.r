# ==============================================================================
# Script Title: 00_get_fastqc.r
# Description: This script gets the quality of the sequences using FastQC
# ==============================================================================
# Author/s: @DiegoFE94 from Github (GitHub)
# Affiliation: Machine Learning for Life Sciences Laboratory (MALL)
# Email: diego.fedreira@udc.es
# Date Created: 2024 (Last Update on January 2025)
# EXTERNAL REQUIREMENT: FASTQC installed and alias in your env/system
# ==============================================================================
# 0. Set Up
# ==============================================================================
# 0.1 Clean Environment
rm(list = ls())
set.seed(1965)

# 0.2 Load Packages
# Uncomment and add your packages here
# library()

# 0.3 Declare or Load Custom Functions
# source("")

# 0.4 Declare Variables
# Inputs
input_path <- "extdata/FASTQ"

# Outputs

# Arguments

# ==============================================================================
# 1.Run FastQC
# ==============================================================================
# 1.1 List the files in the input path
fnfs <- sort(list.files(path, input_path = ".fastq.gz", full.names = TRUE))

# 1.2 make command line for each raw ()
l <- list()
for (i in seq_along(fnfs)){
  l[[i]] <- paste("Fastqc", fnfs[i])
}

# 1.3 Run those commands commmand
for (i in seq_along(fnFs)){
  system(l[[i]])
}

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()