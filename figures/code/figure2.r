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

# Metabolic Syndrome
# A.1
ps <- readRDS("figures/data/phy_genus_MS.rds")
model <- '~Status + GS + BMI'
alpha <- 0.05
linda.res <- linda(phyloseq.obj = ps,
                   feature.dat.type = "count",
                   formula = model,
                   zero.handling = "pseudo-count",
                   alpha = alpha,
                   p.adj.method = "none",
                   prev.filter = 0, 
                   mean.abund.filter = 0)
linda_result <- linda.res$output$StatusMS
df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))
df2plot <- cbind(as(df2plot, "data.frame"),
                as(tax_table(ps)[rownames(df2plot), ], "matrix"))
# Save summaary all genera
write_xlsx(df2plot, path = "figures/results/linda_genus_ms.xlsx")
show_genus <- c("Blautia", "Tyzzerella", "Dorea", "Streptococcus",
                "Family XIII AD3011 group", "Paludicola",
                "[Eubacterium] eligens group", "Prevotellaceae NK3B31 group"
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
                                    size=1))

ggsave(
  p1,
  filename = "fig2a1.svg",
  device = "svg",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")


# Metabolic Syndrome
# A.2
palette <- readRDS(file = "figures/data/custom_palette_species.rds")
ps <-  readRDS("figures/data/phy_species_MS.rds")
model <- '~Status + GS + BMI'
alpha <- 0.05
linda.res <- linda(phyloseq.obj = ps,
                   feature.dat.type = "count",
                   formula = model,
                   zero.handling = "pseudo-count",
                   alpha = alpha,
                   p.adj.method = "none",
                   prev.filter = 0, 
                   mean.abund.filter = 0)
linda_result <- linda.res$output$StatusMS

df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj'))

df2plot <- cbind(as(df2plot, "data.frame"),
                 as(tax_table(ps)[rownames(df2plot), ], "matrix"))
# Save summaary all spps
write_xlsx(df2plot, path = "figures/results/linda_species_ms.xlsx")
# Create a volcano  plot
p2 <- ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
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
                                    size=1))

ggsave(
  p2,
  filename = "fig2b1.svg",
  device = "svg",
  path = "figures/plots/",
  width = 110, 
  height = 69, 
  units = "mm")
