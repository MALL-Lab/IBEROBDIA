## 16S Fastq Processing
library(dada2)
path = "projects/IBEROBDIA/Extdata/FASTQ/"
fnFs <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#plotQualityProfile(fnFs[1])
## Place filtered files in filtered subdirectory
filtFs <- file.path(path, "Filtered_FASTQ", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

#Applying quality filters on sequences. #trimleft= 15 recomended for DADA2 processing ion data
out <- filterAndTrim(fwd = fnFs, filt =filtFs, truncLen= 251,
              maxN=0, maxEE=3, truncQ=2, rm.phix=TRUE, trimLeft=15, 
              compress=TRUE, multithread=TRUE)
head(out)

#Learning errors from sequences 
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+9)
#plotErrors(errF, nominalQ=TRUE)

# Derep Step check (Optional, applied on big data to avoid times computing)
# derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
# names(derepFs) <- sample.names

#Aplying core algorithm to infer real biological sequences. Homology gap and band size parameters recomended for DADA2 processing ion data
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
dadaFs[[1]]
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Removing Quimeras from sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))
#Showing evolution of sequences from raw to final step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
taxa <- assignTaxonomy(seqtab.nochim, "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/tax/silva_nr99_v138.1_train_set.fa.gz", tryRC=TRUE, multithread=TRUE)
taxa <- addSpecies(taxa, "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/tax/silva_species_assignment_v138.1.fa.gz") 

#Save OTU table and Taxa table
saveRDS(object = taxa, file = paste0(path,"tax_table.rds"))
saveRDS(object = seqtab, file = paste0(path,"otu_table.rds"))
