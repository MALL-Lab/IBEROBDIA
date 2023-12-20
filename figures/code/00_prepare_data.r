library(phyloseq)
# Declare variables 
level <- "Species"
experiment <- "MS"
levels <- c("No_MS", "MS")
# Load phyloseq
phy <- readRDS(file = paste0("02_preprocess/data/phy_",
                             level,".rds"))
# Change name of the column experiment as Status to perform the analysis
colnames(phy@sam_data)[colnames(phy@sam_data) == experiment] <- "Status"
phy@sam_data$Status <-factor(x = as.factor(phy@sam_data$Status),
                             levels = levels)
table(phy@sam_data$Status)
# Change genus for "species" or "genus" if change level
saveRDS(object = phy,
        file = paste0("figures/data/phy_species_",experiment,".rds")) 
