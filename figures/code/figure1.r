library(phyloseq)
library(viridis)
library(forcats)
library(dplyr)
library(tibble)
library(ggpubr)
library(ggplot2)
library(phylosmith)
library(writexl)
set.seed(111)
palette <- readRDS(file = "figures/data/custom_palette.rds")

#### A ####
# 1.Phylum Abundancies
phy <- readRDS(file = "02_preprocess/data/phy_Family.rds")
phy <- tax_glom(physeq = phy, taxrank = "Phylum")
phy <- transform_sample_counts(phy, function(ASV) ASV/sum(ASV))
df <- psmelt(phy)
df_top10 <- df %>% 
  group_by(Phylum) %>% 
  summarize(total_abundancia = sum(Abundance)) %>% 
  top_n(10, total_abundancia)
p1 <- df_top10 %>%
  mutate(Phylum = fct_reorder(Phylum, desc(total_abundancia))) %>%
  ggplot(aes(x = Phylum, y = total_abundancia, fill = Phylum)) +
  geom_bar(stat = "identity") +   theme_bw() +
  theme(axis.text.x = element_text(angle = 30,
                                   size = 6,
                                   hjust = 1),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.box = "horizontal",
        legend.text = element_text(size = 7)) +
  scale_y_continuous(limits=c(0, 45)) +
  labs(x = "", y = "Relative Abundance (%)") +
  scale_fill_manual(values = palette)

ggsave(
  p1,
  filename = "fig1a1.svg",
  device = "svg",
  path = "figures/plots/",
  width = 70, 
  height = 69, 
  units = "mm")

# 1.Family Abundancies
phy <- readRDS(file = "02_preprocess/data/phy_Family.rds")
phy <- tax_glom(physeq = phy, taxrank = "Family")
phy <- transform_sample_counts(phy, function(ASV) ASV/sum(ASV))
df <- psmelt(phy)
df_top10 <- df %>% 
  group_by(Family) %>% 
  summarize(total_abundancia = sum(Abundance)) %>% 
  top_n(10, total_abundancia)
df_resto <- df %>% 
  anti_join(df_top10, by = "Family") %>% 
  summarize(Family = "Others", total_abundancia = sum(Abundance))
df_plot <- bind_rows(df_top10, df_resto)
p2 <- df_plot %>%
  mutate(Family = fct_reorder(Family, desc(total_abundancia))) %>%
  ggplot(aes(x = Family, y = total_abundancia, fill = Family)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 35,
                                   size = 6,
                                   hjust = 1),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.box = "horizontal",
        legend.text = element_text(size = 7)) +
  scale_y_continuous(limits=c(0, 45)) +
  labs(x = "", y = "Relative Abundance (%)") +
  scale_fill_manual(values = palette)

ggsave(
  p2,
  filename = "fig1a2.svg",
  device = "svg",
  path = "figures/plots/",
  width = 75, 
  height = 71, 
  units = "mm")

# 1.Genus Abundancies
phy <- readRDS(file = "02_preprocess/data/phy_Genus.rds")
phy <- tax_glom(physeq = phy, taxrank = "Genus")
phy <- transform_sample_counts(phy, function(ASV) ASV/sum(ASV))
df <- psmelt(phy)
df_top10 <- df %>% 
  group_by(Genus) %>% 
  summarize(total_abundancia = sum(Abundance)) %>% 
  top_n(10, total_abundancia)
df_resto <- df %>% 
  anti_join(df_top10, by = "Genus") %>% 
  summarize(Genus = "Others", total_abundancia = sum(Abundance))
df_plot <- bind_rows(df_top10, df_resto)
p3 <- df_plot %>%
  mutate(Genus = fct_reorder(Genus, desc(total_abundancia))) %>%
  ggplot(aes(x = Genus, y = total_abundancia, fill = Genus)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30,
                                   size = 6,
                                   hjust = 1),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.box = "horizontal",
        legend.text = element_text(size = 7)) +
  scale_y_continuous(limits=c(0, 45)) +
  labs(x = "", y = "Relative Abundance (%)") +
  scale_fill_manual(values = palette)

ggsave(
  p3,
  filename = "fig1a3.svg",
  device = "svg",
  path = "figures/plots/",
  width = 70, 
  height = 75, 
  units = "mm")
####A####

####B####
c_palette <- readRDS(file = "figures/data/custom_palette_species.rds")
# B1
phy <- readRDS(file = "02_preprocess/data/phy_Family.rds")
phy <- tax_glom(physeq = phy, taxrank = "Phylum")

sample_data(phy)$MS <- as.factor(sample_data(phy)$MS)
levels(sample_data(phy)$MS)
sample_data(phy)$MS <- relevel(sample_data(phy)$MS, ref = "No_MS")
levels(sample_data(phy)$MS)


target <- "MS"
phy <- merge_samples(x = phy, group = target, fun = sum)
p_ra <- phylogeny_profile(phyloseq_obj = phy, relative_abundance = TRUE,
                         classification = "Phylum", grid = FALSE, merge = TRUE,
                         )
p4 <- p_ra + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   angle = 25,
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.text = element_text(size=6, color = "black"),
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) +  scale_fill_manual(values = c_palette)

exc <-as.data.frame(p4$data[,1:4])
ggsave(
  p4,
  filename = "fig1b1.svg",
  device = "svg",
  path = "figures/plots/",
  width = 100, 
  height = 75, 
  units = "mm")

# B2
phy <- readRDS(file = "02_preprocess/data/phy_Family.rds")
phy <- tax_glom(physeq = phy, taxrank = "Phylum")

sample_data(phy)$DT2_P_H <- as.factor(sample_data(phy)$DT2_P_H)
levels(sample_data(phy)$DT2_P_H)
sample_data(phy)$DT2_P_H <- factor(sample_data(phy)$DT2_P_H, levels = c("Healthy", "PreDT2", "DT2"))
levels(sample_data(phy)$DT2_P_H)

target <- "DT2_P_H"
phy <- merge_samples(x = phy, group = target, fun = sum)
p_ra <- phylogeny_profile(phyloseq_obj = phy, relative_abundance = TRUE,
                          classification = "Phylum", grid = FALSE, merge = TRUE,
)
p5 <- p_ra + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   angle = 25,
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.text = element_text(size=6, color = "black"),
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) +  scale_fill_manual(values = c_palette)
exc <- rbind(exc, as.data.frame(p5$data[,1:4]))

write_xlsx(exc,"figures/results/phylum_ratios.xlsx")

ggsave(
  p5,
  filename = "fig1b2.svg",
  device = "svg",
  path = "figures/plots/",
  width = 100.75, 
  height = 75, 
  units = "mm")

####C####
#### MS ####
ps <- readRDS("figures/data/phy_genus_MS.rds")
colrs <- c("#B3DE69", "#FDB462")
# Alpha Div
rich <- estimate_richness(ps)
df <- data.frame(ps@sam_data)
shanon <- as.numeric(rich$Shannon)
simpson <-as.numeric(rich$Simpson)
status <- as.character(df$Status)
totest <- data.frame(cbind(shanon, simpson, status))
totest$shanon <- as.double(totest$shanon)
totest$simpson <- as.double(totest$simpson)
# Shanon
totest$dummy <- "Shannon"
stat_test <- compare_means(
  shanon ~ status, data = totest,
  method = "wilcox.test",
  p.adjust.method = "fdr",
)
stat_test <- stat_test %>%
  mutate(my_p = paste0("p.adj = ", round(p.adj, digits = 2))) %>%
  mutate(y.position = c(3.5))

p6 <- ggscatter(data = totest, x = "status", y ="shanon",color = "status") +
  geom_boxplot(fill = colrs, color = "azure4",
               outlier.color = "black", outlier.shape = 2 )+
  stat_pvalue_manual(data = stat_test,label = "my_p", size = 2,bracket.size = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   angle = 25,
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + facet_grid(. ~ dummy) + scale_color_manual(values = colrs)

ggsave(
  p6,
  filename = "fig1c1.svg",
  device = "svg",
  path = "figures/plots/",
  width = 42.75, 
  height = 60, 
  units = "mm")

# Simpson
totest$dummy <- "Simpson"
stat_test <- compare_means(
  simpson ~ status, data = totest,
  method = "wilcox.test",
  p.adjust.method = "fdr",
)
stat_test <- stat_test %>%
  mutate(my_p = paste0("p.adj = ", p.adj)) %>%
  mutate(y.position = c(0.95))

p7 <- ggscatter(data = totest, x = "status", y ="simpson",color = "status") +
  geom_boxplot(fill = colrs, color = "azure4",
               outlier.color = "black", outlier.shape = 2 )+
  stat_pvalue_manual(data = stat_test,label = "my_p", size = 2,bracket.size = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   angle = 25,
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + facet_grid(. ~ dummy) + scale_color_manual(values = colrs)

ggsave(
  p7,
  filename = "fig1c2.svg",
  device = "svg",
  path = "figures/plots/",
  width = 42.75, 
  height = 60, 
  units = "mm")

# BetaDiv
# Bray Distance
dist = phyloseq::distance(ps,
                          method ="bray",
                          weighted = F)
ordination = ordinate(ps, method="PCoA", distance=dist)
ps@sam_data$dummy <- "Bray"
set.seed(111)
a <- adonis2(dist ~ sample_data(ps)$Status,permutations = 1000)
a <- a$`Pr(>F)`[1]
p8 <- plot_ordination(ps,type = "samples", ordination, color="Status") +
  geom_point(size = 4.5) +
  annotate(geom="text", x=0, y=0.50,size =2,
           label=paste0("p.adj = ", round(a,digits = 4)),
          color="black") +
  stat_ellipse(aes(fill = Status)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(color = "black",
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = colrs) + facet_grid(. ~ dummy)

ggsave(
  p8,
  filename = "fig1d1.svg",
  device = "svg",
  path = "figures/plots/",
  width = 65, 
  height = 60, 
  units = "mm")

# Chao distance
dist = phyloseq::distance(ps,
                          method ="chao",
                          weighted = F
                          )
ordination = ordinate(ps, method="PCoA", distance=dist)
ps@sam_data$dummy <- "Chao"
a <- adonis2(dist ~ sample_data(ps)$Status,permutations = 1000)
a <- a$`Pr(>F)`[1]
p9 <- plot_ordination(ps,type = "samples", ordination, color="Status") +
  geom_point(size = 4.5) +
  annotate(geom="text", x=0, y=0.30,size =2,
           label=paste0("p.adj = ", round(a,digits = 4)),
           color="black") +
  stat_ellipse(aes(fill = Status)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(color = "black",
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = colrs) + facet_grid(. ~ dummy)

ggsave(
  p9,
  filename = "fig1d2.svg",
  device = "svg",
  path = "figures/plots/",
  width = 65, 
  height = 60, 
  units = "mm")

#### MS ####



#### DT2 ####
ps <- readRDS("figures/data/phy_genus_DT2_P_H.rds")
colrs <- c("#80B1D3", "#8DD3C7", "#BEBADA")
# Alpha Div
rich <- estimate_richness(ps)
df <- data.frame(ps@sam_data)
shanon <- as.numeric(rich$Shannon)
simpson <-as.numeric(rich$Simpson)
status <- as.character(df$Status)
totest <- data.frame(cbind(shanon, simpson, status))
totest$shanon <- as.double(totest$shanon)
totest$simpson <- as.double(totest$simpson)
# Shanon
totest$dummy <- "Shannon"
stat_test <- compare_means(
  shanon ~ status, data = totest,
  method = "kruskal.test",
  p.adjust.method = "fdr",
)
stat_test <- stat_test %>%
  mutate(my_p = paste0("p.adj = ", p.adj)) %>%
  mutate(y.position = c(3.5))

p10 <- ggscatter(data = totest, x = "status", y ="shanon",color = "status") +
  geom_boxplot(fill = colrs, color = "azure4",
               outlier.color = "black", outlier.shape = 2 )+
  annotate(geom="text", x="PreDT2", y=3.50,size =2,
           label=paste0("p.adj = ", round(stat_test$p.adj,digits = 4)),
           color="black") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   angle = 25,
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + facet_grid(. ~ dummy) + scale_color_manual(values = colrs)
ggsave(
  p10,
  filename = "fig1c3.svg",
  device = "svg",
  path = "figures/plots/",
  width = 42.75, 
  height = 60, 
  units = "mm")

# Simpson
totest$dummy <- "Simpson"
stat_test <- compare_means(
  simpson ~ status, data = totest,
  method = "kruskal.test",
  p.adjust.method = "fdr",
)
stat_test <- stat_test %>%
  mutate(my_p = paste0("p.adj = ", p.adj)) %>%
  mutate(y.position = c(0.96))

p11 <- ggscatter(data = totest, x = "status", y ="simpson",color = "status") +
  geom_boxplot(fill = colrs, color = "azure4",
               outlier.color = "black", outlier.shape = 2 )+
  annotate(geom="text", x="PreDT2", y=0.96,size =2,
           label=paste0("p.adj = ", round(stat_test$p.adj,digits = 4)),
           color="black") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   angle = 25,
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + facet_grid(. ~ dummy) + scale_color_manual(values = colrs)

ggsave(
  p11,
  filename = "fig1c4.svg",
  device = "svg",
  path = "figures/plots/",
  width = 42.75, 
  height = 60, 
  units = "mm")

# BetaDiv
# Bray Distance
dist = phyloseq::distance(ps,
                          method ="bray",
                          weighted = F)
ordination = ordinate(ps, method="PCoA", distance=dist)
ps@sam_data$dummy <- "Bray"
set.seed(111)
a <- adonis2(dist ~ sample_data(ps)$Status,permutations = 1000)
a <- a$`Pr(>F)`[1]
p12 <- plot_ordination(ps,type = "samples", ordination, color="Status") +
  geom_point(size = 4.5) +
  annotate(geom="text", x=0, y=0.92,size =2,
           label=paste0("p.adj = ", round(a,digits = 4)),
           color="black") +
  stat_ellipse(aes(fill = Status)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(color = "black",
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = colrs) + facet_grid(. ~ dummy)

ggsave(
  p12,
  filename = "fig1d3.svg",
  device = "svg",
  path = "figures/plots/",
  width = 65, 
  height = 60, 
  units = "mm")

# Chao distance
dist = phyloseq::distance(ps,
                          method ="chao",
                          weighted = F
)
ordination = ordinate(ps, method="PCoA", distance=dist)
ps@sam_data$dummy <- "Chao"
a <- adonis2(dist ~ sample_data(ps)$Status,permutations = 1000)
a <- a$`Pr(>F)`[1]
p13 <- plot_ordination(ps,type = "samples", ordination, color="Status") +
  geom_point(size = 4.5) +
  annotate(geom="text", x=0, y=0.57, size =2,
           label=paste0("p.adj = ", round(a,digits = 4)),
           color="black") +
  stat_ellipse(aes(fill = Status)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(color = "black",
                                   size = 6),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = colrs) + facet_grid(. ~ dummy)

ggsave(
  p13,
  filename = "fig1d4.svg",
  device = "svg",
  path = "figures/plots/",
  width = 65, 
  height = 60, 
  units = "mm")
####DT2####