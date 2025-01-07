# ==============================================================================
# Script Title: 02_linda_analysis.r
# Description: This script performs a LINDA analysis
# ==============================================================================
# Author/s: @DiegoFE94 (GitHub)
# Affiliation: Machine Learning for Life Sciences Laboratory (MALL)
# Email: diego.fedreira@udc.es
# Date Created: 2024 (Last Update on January 2025)
# ==============================================================================
# 0. Set Up
# ==============================================================================
# 0.1 Clean Environment
rm(list = ls())
set.seed(123)

# 0.2 Load Packages
library(phyloseq)
library(MicrobiomeStat)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(microbiome)
library(microViz)
library(ComplexHeatmap)
library(viridisLite)

# 0.3 Declare or Load Custom Functions
# source("")

# 0.4 Declare Variables
# Inputs
input_path <- "02_preprocess/data/"

# Outputs
output_path <- "02_preprocess/data/"
results_path <- "02_preprocess/results/"

# Arguments
level <- "Species"
experiment <- "DT2_P_H"
levels <- c("PreDT2", "Healthy", "DT2")
model <- "~Status + GS + BMI + MS" # Ms: "~Status + GS + BMI" , DT2: "~Status + GS + BMI + MS"
id <- "PreDT2vsDT2"
a <- 0.05

# ==============================================================================
# 1. Load and prepare the data
# ==============================================================================
# 1.1 Load the phyloseq object
phy <- readRDS(file = paste0(input_path, "phy_", level, ".rds"))

# 1.2 Change name of the column experiment as Status to perform the analysis
colnames(phy@sam_data)[colnames(phy@sam_data) == experiment] <- "Status"
phy@sam_data$Status <- factor(x = as.factor(phy@sam_data$Status),
                              levels = levels)
table(phy@sam_data$Status)
target <- "Status"

# ==============================================================================
# 2. Perform the LINDA analysis
# ==============================================================================
# 2.1 Perform the LINDA method
linda.res <- linda(phyloseq.obj = phy,
                   feature.dat.type = "count",
                   formula = model,
                   zero.handling = "pseudo-count",
                   alpha = a,
                   p.adj.method = "none",
                   prev.filter = 0,
                   mean.abund.filter = 0)
linda_result <- linda.res$output$StatusDT2 # We hace to change this for each

# 2.2 Convert results into a dataframe
df2plot <- linda_result %>%
  dplyr::select(c("log2FoldChange", "pvalue", "padj")) %>%
  rownames_to_column("names")

# ==============================================================================
# 3. Visualize results in a volcano plot
# ==============================================================================
ggplot(df2plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$padj > 0.2, "gray",
                                      ifelse(df2plot$pvalue <= 0.05, "red",
                                             ifelse(df2plot$pvalue <= 0.15, "green","green")))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df2plot, pvalue <= 0.05),
                  aes(label = names), size = 2.5, vjust = 1.1, hjust = 1.1) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "-log10(p-value)")  +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    size = 1))

# ==============================================================================
# 4. Procces the results and save them
# ==============================================================================
# 4.1 Select criteria to exclude
sigtab <- linda_result %>%
  dplyr::select(c("log2FoldChange", "pvalue", "padj", "reject")) %>%
  mutate(log2FoldChange = abs(log2FoldChange)) %>%
  filter(pvalue < 0.05)
sigtab <- cbind(as(sigtab, "data.frame"),
                as(tax_table(phy)[rownames(sigtab), ], "matrix"))

# 4.2 Save phyloseq with only taxones of interest
pruned_pseq <- prune_taxa(x = phy, taxa = rownames(sigtab))
pruned_pseq@tax_table
saveRDS(object = pruned_pseq,
        file = paste0(output_path, "pruned_", level, "_", id, "_pseq.rds"))

# 4.3 Save formnated results in csv
sigtab <- sigtab %>%
  dplyr::select(c("log2FoldChange", "pvalue", "padj",
                  "reject", "Family", "Genus", "Species"))
write.csv2(x = sigtab,
           file = paste0(results_path, "results_", level, "_", id, ".csv"),
           sep = ";", col.names = TRUE, row.names = TRUE)

# ==============================================================================
# Session Information
# ==============================================================================
sessionInfo()

# Clean up the environment before restarting the session
rm(list = ls())

# Restart the R session (optional, only if needed)
# .rs.restartR()