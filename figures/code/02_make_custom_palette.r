#Make big Dictionary
#Phylum level
phy <- readRDS(file = "02_preprocess/data/phy_Family.rds")
phy <- tax_glom(physeq = phy, taxrank = "Phylum")
phylums <- as.character(phy@tax_table[,"Phylum"])

# Load data from which extract families and genus
d1 <-  readRDS("figures/data/MS_Genus.rds")
d2 <-  readRDS("figures/data/HP_Genus.rds")
d3 <-  readRDS("figures/data/HD_Genus.rds")
d4 <-  readRDS("figures/data/PD_Genus.rds")
d_all <- rbind(d1,d2,d3,d4)
tf <- readRDS("figures/data/TF.rds")
tg <- readRDS("figures/data/TG.rds")
#Family level
families <- unique(c(d_all$up_level, tf))

#Genus level
genus <- unique(c(d_all$down_level, tg))

# Create colors from existing palettes 
phylum_colors <- brewer.pal(7, "Set3")
# Generate colors for families
family_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(families))
# Generate colors for genus
genus_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(genus))

# Create a named vector using the taxonomic levels and colors
taxonomic_levels <- c(phylums, families, genus)
colors <- c(phylum_colors, family_colors, genus_colors, "#999999")
palette <- setNames(colors, c(taxonomic_levels, "Others"))
saveRDS(object = palette, file ="figures/data/custom_palette.rds")

# Add Species to the palette
palette <- readRDS(file = "figures/data/custom_palette.rds")
d1 <- readRDS("figures/data/MS_Species.rds")
d2 <- readRDS("figures/data/HP_Species.rds")
d3 <- readRDS("figures/data/HD_Species.rds")
d4 <- readRDS("figures/data/PD_Species.rds")
d_all <- rbind(d1,d2,d3,d4)
#Genus level
genus <- unique(c(d_all$up_level))
p <- genus %in% names(palette)
genus <- genus[!p]
genus_colors <- colorRampPalette(brewer.pal(12, "Set2"))(length(genus))
#Species level
species <- unique(c(d_all$down_level))
species_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(species))
#Create vector palete
taxonomic_levels <- c(genus, species)
colors <- c(genus_colors, species_colors)
palette2 <- setNames(colors, taxonomic_levels)
palette_species <- c(palette, palette2)
saveRDS(object = palette_species, file ="figures/data/custom_palette_species.rds")
