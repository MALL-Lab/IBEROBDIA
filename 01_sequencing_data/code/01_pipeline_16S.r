# DADA2 16S pipeline
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/IBEROBDIA')
source(file = "01_sequencing_data/code/config_file.r")

## 16S Fastq Processing
library(dada2)
fnFs <- sort(list.files(input.dir.path, pattern=".fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#plotQualityProfile(fnFs[1])

## Place filtered files in filtered subdirectory
filter.path = paste(input.dir.path, "Filtered_FASTQ", sep= "/")
if (dir.exists(filter.path) == FALSE) {
  dir.create(filter.path)
  message(paste('Creating ', filter.path, ' directory!'))
}
filtFs <- file.path(filter.path, paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

#Applying quality filters on sequences. #trimleft= 15 recomended for DADA2 processing ion data
out <- filterAndTrim(fwd = fnFs, filt = filtFs, truncLen = dada2.truncLen, maxN = dada2.maxN,
                     maxEE = dada2.maxEE, truncQ = dada2.truncQ, trimLeft = dada2.trimLeft,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

#Learning errors from sequences 
errF <- learnErrors(filtFs, multithread=TRUE, nbases = dada2.nbases)
#plotErrors(errF, nominalQ=TRUE)

# Derep Step check (Optional, applied on big data to avoid times computing)
# derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
# names(derepFs) <- sample.names

#Aplying core algorithm to infer real biological sequences. Homology gap and band size parameters recomended for DADA2 processing ion data
dadaFs <- dada(filtFs, err = errF, multithread = TRUE,
               HOMOPOLYMER_GAP_PENALTY = dada2.homopolymer_gap_penalty, BAND_SIZE = dada2.band_size)
dadaFs[[1]]
seqtab <- makeSequenceTable(dadaFs)
print("Dimensions of sequence table before quimeras removal:")
dim(seqtab)
d1 = dim(seqtab)
# Inspect distribution of sequence lengths
print("Distribution of sequence lengths before quimeras removal:")
table(nchar(getSequences(seqtab)))

#Removing Quimeras from sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method = dada2.method, multithread = TRUE, verbose = TRUE)
print("Dimensions of sequence table after quimeras removal:")
dim(seqtab.nochim)
d2 = dim(seqtab.nochim)
print("Distribution of sequence lengths after quimeras removal:")
table(nchar(getSequences(seqtab.nochim)))
k = sum(seqtab.nochim)/sum(seqtab)
d3 = d2[2]/d1[2]
print(paste((1 -round(d3,digits = 3))*100,"% of the sequences were quimeras" ))
print(paste("The abundance of these quimeras only represent", round(x = k, digits = 3)*100, "% of the total abundance"))

#Showing evolution of sequences from raw to final step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
print("Summary of the first 5 samples along the Pipeline:")
head(track)

# Taxonomical assignment
taxa <- assignTaxonomy(seqs = seqtab.nochim, refFasta = dada2.path.refFasta,
                       tryRC = dada2.tryRC, multithread = TRUE)
taxa <- addSpecies(taxtab = taxa, refFasta = dada2.path.refFasta.species,
                   tryRC =  dada2.tryRC.species) 

#Save OTU table and Taxa table
saveRDS(object = taxa, file = paste0(out.path,"tax_table.rds"))
saveRDS(object = seqtab, file = paste0(out.path,"otu_table.rds"))
