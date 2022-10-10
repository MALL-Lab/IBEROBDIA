library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Load Sequences from DADA2
path ="../data/"
files = list.files(path = path, pattern = "TL")
files = files[grep('251', files)]

taxas = list()
trees = list()
for (i in seq_along(files)){
  # Load tax table
  taxas[[i]] = readRDS(file = paste0(path,"/", files[i],"/tax_table.rds" ))
  
  # Extract sequences from DADA2 output
  sequences<-getSequences(taxas[[1]])
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
  
}







