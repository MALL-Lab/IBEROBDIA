# ==============================================================================
# Script Title: 01_pipeline_16S.r
# Description: This script is used to process 16S rRNA gene sequencing data.
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
# 0.1.1 Script launched in a HPC cluster
#setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/IBEROBDIA')

# 0.2 Load Packages
# Uncomment and add your packages here
library(dada2)

# 0.3 Declare or Load Custom Functions
getN <- function(x) sum(getUniques(x))

# 0.3.1 Load the configuration file
source(file = "01_sequencing_data/code/config_file.r")

# 0.4 Declare Variables

# Inputs

# Outputs
filter.path = paste(input.dir.path, "Filtered_FASTQ", sep= "/")

# Arguments

# ==============================================================================
# 1. Prepare data
# ==============================================================================
# 1.1 List the files in the input path
fnFs <- sort(list.files(input.dir.path, pattern = ".fastq.gz",
                        full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# 1.2 Create the output directory
if (dir.exists(filter.path) == FALSE) {
  dir.create(filter.path)
  message(paste("Creating ", filter.path, " directory!"))
}
filtFs <- file.path(filter.path, paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

# ==============================================================================
# Applying quality filters on sequences.
# ==============================================================================
out <- filterAndTrim(fwd = fnFs, filt = filtFs, truncLen = dada2.truncLen,
                     maxN = dada2.maxN, maxEE = dada2.maxEE,
                     truncQ = dada2.truncQ, trimLeft = dada2.trimLeft,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

# ==============================================================================
# 3. Learning errors from sequences
# ==============================================================================
errF <- learnErrors(filtFs, multithread = TRUE, nbases = dada2.nbases)

# ==============================================================================
# 4. Aplying core algorithm to infer real biological sequences.
# ==============================================================================
# 4.1 Apply DADA2 to the filtered sequences
dadaFs <- dada(filtFs, err = errF, multithread = TRUE,
               HOMOPOLYMER_GAP_PENALTY = dada2.homopolymer_gap_penalty,
               BAND_SIZE = dada2.band_size)
dadaFs[[1]]
# 4.2 Construct the sequence table
seqtab <- makeSequenceTable(dadaFs)
print("Dimensions of sequence table before quimeras removal:")
dim(seqtab)
d1 <- dim(seqtab)

# 4.3 Inspect distribution of sequence lengths
print("Distribution of sequence lengths before quimeras removal:")
table(nchar(getSequences(seqtab)))

# ==============================================================================
# 5. Removing Quimeras from sequences
# ==============================================================================
# 5.1 Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = dada2.method, multithread = TRUE, verbose = TRUE)

# 5.2 Inspect the results
print("Dimensions of sequence table after quimeras removal:")
dim(seqtab.nochim)
d2 <- dim(seqtab.nochim)
print("Distribution of sequence lengths after quimeras removal:")
table(nchar(getSequences(seqtab.nochim)))
k <- sum(seqtab.nochim) / sum(seqtab)
d3 <- d2[2] / d1[2]
print(paste((1 - round(d3, digits = 3)) * 100,
            "% of the sequences were quimeras"))
print(paste("The abundance of these quimeras only represent",
            (1 - round(x = k, digits = 3)) * 100, "% of the total abundance"))

# ==============================================================================
# 6. Showing evolution of sequences from raw to final step
# ==============================================================================
# 6.1 Create a track table
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
print("Summary of the first 5 samples along the Pipeline:")
head(track)

# ==============================================================================
# 7. Taxonomical assignment
# ==============================================================================
# 7.1 Assign taxonomy
taxa <- assignTaxonomy(seqs = seqtab.nochim, refFasta = dada2.path.refFasta,
                       tryRC = dada2.tryRC, multithread = TRUE)

# 7.2 Add species to the taxa table
taxa <- addSpecies(taxtab = taxa, refFasta = dada2.path.refFasta.species,
                   tryRC =  dada2.tryRC.species)

# ==============================================================================
# 8. Save ASV and Tax tables
# ==============================================================================
saveRDS(object = taxa, file = paste0(out.path, "tax_table.rds"))
saveRDS(object = seqtab.nochim, file = paste0(out.path, "otu_table.rds"))

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()