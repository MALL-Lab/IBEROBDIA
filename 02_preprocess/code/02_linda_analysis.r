# Script to make Linda analysis 
# -------------------------------------
rm(list = ls())
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

set.seed(123)
make_circos_plot <- function(df, group, col.pal){
  require(circlize)
  require(RColorBrewer)
  require(grDevices)
  
  circos.clear()
  circos.par(start.degree=90)
  circos.par("canvas.xlim" = c(-2.5, 2.5))  # Ajustar el espacio en blanco horizontal
  circos.par("canvas.ylim" = c(-2.5, 2.5))
  chordDiagram(df, grid.col = col.pal,
               annotationTrack = c('grid'),
               annotationTrackHeight = c(0.03, 0.01),
               group = group, big.gap = 10, small.gap = 0.5
               #col = colorRampPalette(brewer.pal(8, "Set2"))(20)
  )
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    if(sector.name %in% c(unique(df$up_level),unique(df$down_level))){
      circos.text(CELL_META$xcenter,
                  ylim[1] + cm_h(2),
                  sector.name,
                  facing = "clockwise",
                  niceFacing = TRUE,
                  adj = c(0, 0.5),
                  cex = 0.7,
                  #col=grid.col[sector.name],
                  #font = 2
      )
    }
  }, bg.border = NA)
}
# Declare variables 
level <- "Species"
experiment <- "DT2_P_H"
levels <- c("PreDT2", "Healthy", "DT2")
model <- '~Status + GS + BMI + MS' # Ms: '~Status + GS + BMI' , DT2: '~Status + GS + BMI + MS'
id <- "PreDT2vsDT2"
a <- 0.05
# Load phyloseq
phy <- readRDS(file = paste0("02_preprocess/data/phy_",
                             level,".rds"))

# Change name of the column experiment as Status to perform the analysis
colnames(phy@sam_data)[colnames(phy@sam_data) == experiment] <- "Status"
phy@sam_data$Status <-factor(x = as.factor(phy@sam_data$Status),
                             levels = levels)
table(phy@sam_data$Status)
target <- "Status"

linda.res <- linda(phyloseq.obj = phy,
  feature.dat.type = "count",
  formula = model,
  zero.handling = "pseudo-count",
  alpha = a,
  p.adj.method = "none",
  prev.filter = 0, 
  mean.abund.filter = 0)
linda_result <- linda.res$output$StatusDT2 # We hace to change this for each status

df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj')) %>%
  rownames_to_column("names")

# Create a volcano  plot
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
                                    fill=NA,
                                    size=1))

# select criteri to exclude
sigtab <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj','reject')) %>%
  mutate(log2FoldChange = abs(log2FoldChange)) %>%
  filter(pvalue < 0.05)

sigtab <- cbind(as(sigtab, "data.frame"),
                as(tax_table(phy)[rownames(sigtab), ], "matrix"))

# Save phyloseq with only taxones of interest
pruned_pseq <- prune_taxa(x = phy, taxa = rownames(sigtab))
pruned_pseq@tax_table
saveRDS(object = pruned_pseq,
        file = paste0("02_preprocess/data/pruned_", level, "_",id, "_pseq.rds"
                      ))
# Save formnated results in csv
sigtab <- sigtab %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj',
                  'reject', 'Family', 'Genus', 'Species'))
write.csv2(x = sigtab,file = paste0("02_preprocess/results/results_",level,"_",
                                    id,".csv"),
           sep = ";", col.names = TRUE,row.names = TRUE)

# Circos plot
## Prepare data
if (level == "Genus") {
  circus <- sigtab %>%
    dplyr::select(Family, Genus, log2FoldChange) %>%
    dplyr::rename(up_level = Family, down_level = Genus) %>%
    arrange(down_level)
  rownames(circus)<-NULL
  nm <- unique(c(circus$down_level, circus$up_level))
  my_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(nm))
  grid.col <- data.frame(nm, my_colors)
  my_palette <- tibble::deframe(grid.col)
  links <- structure(c(rep('down_level',length(unique(circus$down_level))),
                       rep('up_level',length(nm)-length(unique(circus$down_level)))),
                     names = nm)
} else if (level == "Species") {
  circus <- sigtab %>% 
    select(Genus, Species, log2FoldChange) %>%
    dplyr::rename(up_level = Genus, down_level = Species) %>%
    arrange(down_level)
  rownames(circus)<-NULL
  nm <- unique(c(circus$down_level, circus$up_level))
  my_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(nm))
  grid.col <- data.frame(nm, my_colors)
  my_palette <- tibble::deframe(grid.col)
  links <- structure(c(rep('down_level', length(unique(circus$down_level))),
                       rep('up_level',length(nm)-length(unique(circus$down_level)))),
                     names = nm)
} else {
  stop("Please insert a valid level (Species or Genus)")
}
## Print plot
saveRDS(object = circus, file = "figures/data/PD_Species.rds") # Save name to make a palette
make_circos_plot(df = circus, group = links , col.pal = my_palette)
