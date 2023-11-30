# Load packages
library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)
TL <- "TL_251"
path.phy = paste0("01_sequencing_data/data/",TL,"/",TL,"_phy_v2.rds")
# Load phy.object
phy = readRDS(path.phy)
rank = 5
# Summary
ntaxa(phy)
nsamples(phy)
sample_names(phy)[1:5]  
rank_names(phy)  
sample_variables(phy)  
otu_table(phy)[1:5, 1:5]
tax_table(phy)[1:5, 1:7]

#1. Agglomerate in genus (There are very few species, if we agglomerate by these we lose most of them)
g_phy = tax_glom(physeq = phy,taxrank=rank_names(phy)[rank])

#2. Filter out ASVs that do not meet certain criteria.
# Remove taxa that are not seen more than 3 times in at least 20% of the samples.
gp = filter_taxa(g_phy, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# Remove samples that do not have a minimum of 20 counts. - In this case al lsamples pass the criteria- 
gp = prune_samples(sample_sums(gp)>=20, gp)
any(taxa_sums(phy) == 0)
saveRDS(object = gp, file = paste0("02_preprocess/data/phy_",rank_names(phy)[rank],"_v2.rds"))
