library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
# Load Sequences from DADA2
taxa =readRDS("projects/IBEROBDIA/Data/TL_251/tax_table.rds")

# Extract sequences from DADA2 output
sequences<-getSequences(taxa)
names(sequences)<-sequences

# Run Sequence Alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

# Change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

# Create distance matrix
dm <- dist.ml(phang.align)

# Perform Neighbor joining
treeNJ <- NJ(dm) # Note, tip order != sequence order

# Internal maximum likelihood
fit = pml(treeNJ, data=phang.align)

# Negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR$tree

