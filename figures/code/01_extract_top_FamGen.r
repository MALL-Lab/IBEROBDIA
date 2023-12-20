# Family
phy <- readRDS(file = "02_preprocess/data/phy_Family.rds")
phy <- tax_glom(physeq = phy, taxrank = "Family")
phy <- transform_sample_counts(phy, function(ASV) ASV/sum(ASV))
df <- psmelt(phy)
df_top10 <- df %>% 
  group_by(Family) %>% 
  summarize(total_abundancia = sum(Abundance)) %>% 
  top_n(10, total_abundancia)

tf <- df_top10$Family
saveRDS(object = tf, "figures/data/TF.rds")

# Genus 
phy <- readRDS(file = "02_preprocess/data/phy_Genus.rds")
phy <- tax_glom(physeq = phy, taxrank = "Genus")
phy <- transform_sample_counts(phy, function(ASV) ASV/sum(ASV))
df <- psmelt(phy)
df_top10 <- df %>% 
  group_by(Genus) %>% 
  summarize(total_abundancia = sum(Abundance)) %>% 
  top_n(10, total_abundancia)

tg <- df_top10$Genus
saveRDS(object = tg, "figures/data/TG.rds")