# Abundances
make_clr <- function(phylo){
  # Taxa have to be rows
  library(mixOmics)
  library(phyloseq)
  pseq <- phylo
  otu <- as.data.frame(otu_table(pseq))
  otu <- t(otu)
  otu <- otu + 1e-05
  otu <- logratio.transfo(otu, logratio = 'CLR')
  otu <- t(otu)
  class(otu) <- 'matrix'
  otu_table(pseq) <- otu_table(otu, taxa_are_rows = TRUE)
  return(pseq)
}


library(ggpubr)
palette <- readRDS(file = "figures/data/custom_palette.rds")
ps <- readRDS(file = "02_preprocess/data/phy_Phylum.rds")
taxa_are_rows(ps)
ps_phylum <- phyloseq::tax_glom(ps, "Phylum")
taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
otu_table(ps_phylum)[1:5, 1:5]
#ps_phylum <- make_clr(ps_phylum)
#ps_phylum <- transform_sample_counts(ps_phylum, function(ASV) ASV/sum(ASV))


sample_data(ps_phylum)$MS <- as.factor(sample_data(ps_phylum)$MS)
levels(sample_data(ps_phylum)$MS)
sample_data(ps_phylum)$MS <- relevel(sample_data(ps_phylum)$MS, ref = "No_MS")
levels(sample_data(ps_phylum)$MS)


sample_data(ps_phylum)$DT2_P_H <- as.factor(sample_data(ps_phylum)$DT2_P_H)
levels(sample_data(ps_phylum)$DT2_P_H)
sample_data(ps_phylum)$DT2_P_H <- factor(sample_data(ps_phylum)$DT2_P_H, levels = c("Healthy", "PreDT2", "DT2"))
levels(sample_data(ps_phylum)$DT2_P_H)
levels(sample_data(ps_phylum)$DT2_P_H) <- c("Healthy", "PreT2D", "T2D")




phyloseq::psmelt(ps_phylum) %>%
ggplot(data = ., aes(x = MS, y = Abundance)) +
  geom_boxplot(outliers = FALSE, aes(fill= OTU, alpha= 0.8)) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  stat_compare_means(#method = "t.test",
   label = "p.signif") + theme(legend.position = "none")


my_comparisons <- list( c("Healthy", "PreT2D"),
 c("Healthy", "T2D"),
 c("PreT2D", "T2D")
 )

phyloseq::psmelt(ps_phylum) %>%
ggplot(data = ., aes(x = DT2_P_H, y = Abundance)) +
  geom_boxplot(outliers = FALSE, aes(fill= OTU, alpha= 0.8)) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  stat_compare_means(#method = "t.test",
  comparisons = my_comparisons,
   label = "p.signif") + theme(legend.position = "none")



