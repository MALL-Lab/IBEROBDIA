# Load packages
library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
TL = "TL_251"
path.phy = paste0("../../01_sequencing_data/data/",TL,"/",TL,"_phy.rds")
# Load phy.object
phy = readRDS(path.phy)

# Summary
ntaxa(phy)
nsamples(phy)
sample_names(phy)[1:5]  
rank_names(phy)  
sample_variables(phy)  
otu_table(phy)[1:5, 1:5]  
tax_table(phy)[1:5, 1:7]

p = data.frame(sample_data(phy))
# R adding a column to dataframe based on values in other columns(Health vs Affected(DT2,TGA,GAA):
p<- p %>% 
  mutate(Statusv2 = if_else(Status == "Healthy", "Healthy", "Affected"))
sample_data(phy) <- p
p<- p %>% 
  mutate(Statusv3 = if_else(Status == "DT2", "DT2", "No_DT2"))
sample_data(phy) <- p

# Preprocess Phyloseq
i = 7
#1. Agglomerate in genus (There are very few species, if we agglomerate by these we lose most of them)
g_phy = tax_glom(physeq = phy,taxrank=rank_names(phy)[i])

#2. Filter out ASVs that do not meet certain criteria.
# Remove taxa that are not seen more than 3 times in at least 20% of the samples.
gp = filter_taxa(g_phy, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# Remove samples that do not have a minimum of 20 counts. - In this case al lsamples pass the criteria- 
gp = prune_samples(sample_sums(gp)>=20, gp)

saveRDS(object = gp, file = paste0("../data/phy_",rank_names(phy)[i],".rds"))






