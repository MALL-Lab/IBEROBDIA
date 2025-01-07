library(phyloseq)
library(MicrobiomeStat)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(microViz)
library(ComplexHeatmap)
library(viridisLite)
library("writexl")
palette <- readRDS(file = "figures/data/custom_palette.rds")
colrs <- c("#6194b9", "#55ad9f", "#7770a5")
# Healthy vs Prediabetic
# A1
ps = readRDS("figures/data/phy_genus_DT2_P_H.rds")
model <- '~Status + GS + BMI + MS'
alpha <- 0.05
linda.res <- linda(phyloseq.obj = ps,
                   feature.dat.type = "count",
                   formula = model,
                   zero.handling = "pseudo-count",
                   alpha = alpha,
                   p.adj.method = "none",
                   prev.filter = 0, 
                   mean.abund.filter = 0)
linda_result <- linda.res$output$StatusPreDT2
df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))

df2plot <- cbind(as(df2plot, "data.frame"),
                 as(tax_table(ps)[rownames(df2plot), ], "matrix"))

write_xlsx(df2plot, path = "figures/results/linda_genus_HP.xlsx")

show_genus <- c("Subdoligranulum", "Dorea", "Anaerostipes", "Fusicatenibacter",
                "[Eubacterium] hallii group", "Ligilactobacillus",
                "Coprococcus", "Acidaminococcus", "GCA-900066575", "Collinsella",
                "Bifidobacterium", "Defluviitaleaceae UCG-011", "Faecalibacterium"
                
)
df2plot$show <- df2plot$Genus %in% show_genus


# Create a volcano  plot
p1 <- ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$pvalue > 0.2, "gray",
                                      ifelse(df2plot$pvalue <= 0.05, "red","green"))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df2plot, show == TRUE),
                  aes(label = Genus), size = 2.5, vjust = 1.1, hjust = 1.1) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)")  +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    size=1)) +
  annotate("text", x = -9, y = 0.1, label = "Healthy",
           color = "#6194b9", size = 3, hjust = 0) +
  annotate("text", x = 7.5, y = 0.1, label = "PreT2D",
           color = "#55ad9f", size = 3, hjust = 1)
p1
ggsave(
  p1,
  filename = "fig3a1.png",
  device = "png",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")

# Healthy vs DT2
# B1
linda_result <- linda.res$output$StatusDT2

df2plot <- linda_result %>%
dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))

df2plot <- cbind(as(df2plot, "data.frame"),
                 as(tax_table(ps)[rownames(df2plot), ], "matrix"))

write_xlsx(df2plot, path = "figures/results/linda_genus_HD.xlsx")
show_genus <- c("Subdoligranulum", "Dorea", "Anaerostipes", "Bifidobacterium",
                "[Eubacterium] ventriosum group", "Ligilactobacillus",
                "Acidaminococcus", "GCA-900066575", "Collinsella",
                "Bifidobacterium", "Defluviitaleaceae UCG-011", "Faecalibacterium"
                
)
df2plot$show <- df2plot$Genus %in% show_genus

# Create a volcano  plot
p2 <- ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$padj > 0.2, "gray",
                                      ifelse(df2plot$pvalue <= 0.05, "red","green"))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df2plot, show == TRUE),
                  aes(label = Genus), size = 2.5, vjust = 1.1, hjust = 1.1) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)")  +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    size=1)) + 
  annotate("text", x = -9, y = 0.1, label = "Healthy",
           color = "#6194b9", size = 3, hjust = 0) +
  annotate("text", x = 5.5, y = 0.1, label = "T2D",
           color = "#7770a5", size = 3, hjust = 1)

ggsave(
  p2,
  filename = "fig3b1.png",
  device = "png",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")

# Prediabetic vs DT2
# C1
ps = readRDS("figures/data/phy_genus_DT2_P_H.rds")
levels <- c("PreDT2", "Healthy", "DT2")
ps@sam_data$Status <-factor(x = as.factor(ps@sam_data$Status),
                             levels = levels)
model <- '~Status + GS + BMI + MS'
linda.res <- linda(phyloseq.obj = ps,
                   feature.dat.type = "count",
                   formula = model,
                   zero.handling = "pseudo-count",
                   alpha = alpha,
                   p.adj.method = "none",
                   prev.filter = 0, 
                   mean.abund.filter = 0)
linda_result <- linda.res$output$StatusDT2
df2plot <- 
  df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))
df2plot <- cbind(as(df2plot, "data.frame"),
                 as(tax_table(ps)[rownames(df2plot), ], "matrix"))
write_xlsx(df2plot, path = "figures/results/linda_genus_PD.xlsx")
show_genus <- c("Dialister", "Dorea", "Sutterella", "Prevotella",
                "[Eubacterium] ventriosum group", "[Eubacterium] hallii group",
                "Erysipelotrichaceae UCG-003", "Fusicatenibacter", "Roseburia",
                "CAG-56", "Lachnospiraceae FCS020 group"
                
)
df2plot$show <- df2plot$Genus %in% show_genus

# Create a volcano  plot
p3 <- ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$padj > 0.2, "gray",
                                      ifelse(df2plot$pvalue <= 0.05, "red","green"))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df2plot, show == TRUE),
                  aes(label = Genus), size = 2.5, vjust = 1.1, hjust = 1.1) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)")  +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    size=1)) +
  annotate("text", x = -2.7, y = 0.1, label = "PreT2D",
           color = "#55ad9f", size = 3, hjust = 1) +
  annotate("text", x = 5.5, y = 0.1, label = "T2D",
           color = "#7770a5", size = 3, hjust = 1)

ggsave(
  p3,
  filename = "fig3c1.png",
  device = "png",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")

level = "Species"
c_palette <- readRDS(file = "figures/data/custom_palette_species.rds")


# Healthy vs Prediabetic Species
# A2
ps <- readRDS("figures/data/phy_species_DT2_P_H.rds")
model <- '~Status + GS + BMI + MS'
alpha <- 0.05
linda.res <- linda(phyloseq.obj = ps,
                   feature.dat.type = "count",
                   formula = model,
                   zero.handling = "pseudo-count",
                   alpha = alpha,
                   p.adj.method = "none",
                   prev.filter = 0, 
                   mean.abund.filter = 0)
linda_result <- linda.res$output$StatusPreDT2

df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))

df2plot <- cbind(as(df2plot, "data.frame"),
                 as(tax_table(ps)[rownames(df2plot), ], "matrix"))

write_xlsx(df2plot, path = "figures/results/linda_species_HP.xlsx")

# Create a volcano  plot
p4 <- ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$padj > 0.2, "gray",
                                      ifelse(df2plot$pvalue <= 0.05, "red",
                                             ifelse(df2plot$pvalue <= 0.15, "green","green")))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df2plot, pvalue <= 0.05),
                  aes(label = Species), size = 2.5, vjust = 1.1, hjust = 1.1) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)")  +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    size=1)) +
  annotate("text", x = -5.9, y = 0.1, label = "Healthy",
           color = "#6194b9", size = 3, hjust = 0) +
  annotate("text", x = 5.9, y = 0.1, label = "PreT2D",
           color = "#55ad9f", size = 3, hjust = 1)

ggsave(
  p4,
  filename = "fig3a2.png",
  device = "png",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")

# Healthy vs DT2 Species
linda_result <- linda.res$output$StatusDT2

df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))

df2plot <- cbind(as(df2plot, "data.frame"),
                 as(tax_table(ps)[rownames(df2plot), ], "matrix"))

write_xlsx(df2plot, path = "figures/results/linda_species_HD.xlsx")

# Create a volcano  plot
p5 <- ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$padj > 0.2, "gray",
                                      ifelse(df2plot$pvalue <= 0.05, "red",
                                             ifelse(df2plot$pvalue <= 0.15, "green","green")))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df2plot, pvalue <= 0.05),
                  aes(label = Species), size = 2.5, vjust = 1.1, hjust = 1.1) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)")  +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    size=1)) + 
  annotate("text", x = -5.9, y = 0.1, label = "Healthy",
           color = "#6194b9", size = 3, hjust = 0) +
  annotate("text", x = 3.5, y = 0.1, label = "T2D",
           color = "#7770a5", size = 3, hjust = 1)

ggsave(
  p5,
  filename = "fig3b2.png",
  device = "png",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")

# Prediabetic vs DT2 Species
# C3
ps = readRDS("figures/data/phy_species_DT2_P_H.rds")
levels <- c("PreDT2", "Healthy", "DT2")
ps@sam_data$Status <-factor(x = as.factor(ps@sam_data$Status),
                            levels = levels)
model <- '~Status + GS + BMI + MS'
linda.res <- linda(phyloseq.obj = ps,
                   feature.dat.type = "count",
                   formula = model,
                   zero.handling = "pseudo-count",
                   alpha = alpha,
                   p.adj.method = "none",
                   prev.filter = 0, 
                   mean.abund.filter = 0)
linda_result <- linda.res$output$StatusDT2

df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))

df2plot <- cbind(as(df2plot, "data.frame"),
                 as(tax_table(ps)[rownames(df2plot), ], "matrix"))

write_xlsx(df2plot, path = "figures/results/linda_species_PD.xlsx")

# Create a volcano  plot
p6 <- ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$padj > 0.2, "gray",
                                      ifelse(df2plot$pvalue <= 0.05, "red",
                                             ifelse(df2plot$pvalue <= 0.15, "green","green")))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df2plot, pvalue <= 0.05),
                  aes(label = Species), size = 2.5, vjust = 1.1, hjust = 1.1) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)")  +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    size=1)) +
  annotate("text", x = -2.1, y = 0.1, label = "PreT2D",
           color = "#55ad9f", size = 3, hjust = 1) +
  annotate("text", x = 2.8, y = 0.1, label = "T2D",
           color = "#7770a5", size = 3, hjust = 1)

ggsave(
  p6,
  filename = "fig3c2.png",
  device = "png",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")
