# Script to make DEA analysis 
# -------------------------------------
# 3 variables:
# level :"Genus" or "Species"
#
# experiment:"DT2_All" "DT2_NoDT2 "DT2_P_H" "DT2_H" "Aff_H" "MS" "GS" "BMI" "IR"
# "DPH_MS" "DPH_GS" "DPH_BMI" "DPH_GS_BMI" "AffH_MS" "AffH_GS" "AffH_BMI"
# "AffH_GS_BMI" "AffH_MS_GS_BMI" "MS_GS" "MS_BMI" "MS_IR" "MS_GS_BMI" "BMI_IR"
#
# filt_by:"padj" "pvalue" "logFC" "all"
# -------------------------------------
rm(list = ls())
# Declare variables 
level <- "Genus"
experiment <- "AffH_BMI"
filt_by <- "all"

# Required Packages 
require(phyloseq)
require(DESeq2)
require(ggpubr)
require(tidyverse)
require(tibble)
require(dplyr)
require(RColorBrewer)

# Define functions
# Circos Plot
make_circos_plot <- function(df, group, col.pal){
  library(circlize)
  library(RColorBrewer)
  library(grDevices)
  
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
# Plot DEA (gonna vary depending on the level)
if (level == "Genus") {
  plot_dea <- function(sigtab) {
    require(viridis)
    require(ggplot2)
    # Phylum order
    x <- tapply(sigtab$log2FoldChange,
                sigtab$Family,
                function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Family <- factor(as.character(sigtab$Family), levels = names(x))
    
    # Genus order
    x <- tapply(sigtab$log2FoldChange,
                sigtab$Genus,
                function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Genus <- factor(as.character(sigtab$Genus), levels = names(x))
    
    ggplot(sigtab,
           aes(x = Genus,
               y = log2FoldChange,
               color = Family)) +scale_color_manual(values = viridis(length(unique(sigtab$Family))))+
      geom_point(size = 6) +
      theme_bw() +
      theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = -30,
                                                                      hjust = 0,
                                                                      vjust = 0.5))
    
  }
} else if (level == "Species") {
  plot_dea <- function(sigtab) {
    require(viridis)
    require(ggplot2)
    # Genus order
    x <- tapply(sigtab$log2FoldChange,
                sigtab$Family,
                function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Family <- factor(as.character(sigtab$Family), levels = names(x))
    
    # Species order
    x <- tapply(sigtab$log2FoldChange,
                sigtab$Species,
                function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Species <- factor(as.character(sigtab$Species), levels = names(x))
    
    ggplot(sigtab,
           aes(x = Species,
               y = log2FoldChange,
               color = Family)) +scale_color_manual(values = viridis(length(unique(sigtab$Family))))+
      geom_point(size = 6) +
      theme_bw() +
      theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = -30,
                                                                      hjust = 0,
                                                                      vjust = 0.5))
    
  }
} else {
  stop("Please insert a valid level (Species or Genus)")
}

# Load phyloseq SINO CARGAR EL OTRO Y CAMBIAR LA METADATA
phy <- readRDS(file = paste0("02_preprocess/data/phy_",
                             level,"_v2.rds"))
# Filter some samples in combination of columns taht have 1 or 2 samples or
# samples considered like outliers
if (experiment == "DT2_H") {
  phy <- phy %>%
    subset_samples(DT2_P_H != "PreDT2")
  colnames(phy@sam_data)[colnames(phy@sam_data) == "DT2_P_H"] <- "DT2_H"
} else if (experiment == "DT2_NoDT2"){
  phy@sam_data$DT2_P_H <- ifelse(phy@sam_data$DT2_P_H != "DT2", "No_DT2",
                                 phy@sam_data$DT2_P_H)
  colnames(phy@sam_data)[colnames(phy@sam_data) == "DT2_P_H"] <- "DT2_NoDT2"
} else if (experiment == "DPH_GS"|| experiment == "DPH_GS_BMI"){
  phy <- phy %>%
    subset_samples(DPH_GS != "Healthy_AG")
} else if (experiment == "AffH_GS" || experiment == "AffH_GS_BMI"){
  phy <- phy %>%
    subset_samples(AffH_GS != "Healthy_AG")
} else if (experiment == "AffH_MS_GS_BMI"){
  phy <- phy %>%
    subset_samples(AffH_MS_GS_BMI != "Healthy_MS_AG_OB") %>%
    subset_samples(AffH_MS_GS_BMI != "Healthy_MS_NG_NW")
} else if (experiment == "MS_BMI"|| experiment == "MS_GS_BMI"){
  phy <- phy %>%
    subset_samples(MS_BMI != "MS_NW")
} else if (experiment == "IR" || experiment == "MS_IR"){
  phy <- phy %>%
    subset_samples(IR != "NA")
} else if (experiment == "BMI_IR"){
  phy <- phy %>%
    subset_samples(IR != "NA") %>%
    subset_samples(BMI_IR != "NW_ER") %>%
    subset_samples(BMI_IR != "NW_SR")
}
# Chnage name of the column experiment as Status to perform the analysis
colnames(phy@sam_data)[colnames(phy@sam_data) == experiment] <- "Status"
target <- "Status"

table(phy@sam_data$Status)
# Deseq Analysis
phy@sam_data$Status <- as.factor(phy@sam_data$Status) 
diaggen <- phyloseq_to_deseq2(phy, ~Status)
diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
diaggen <- DESeq(diaggen, test= "Wald", fitType = "parametric")
res <- DESeq2::results(diaggen, cooksCutoff = FALSE)
if (filt_by == "pvalue"){
  sigtab = res[which(res$pvalue < 0.05), ]
} else if (filt_by == "padj"){
  sigtab = res[which(res$padj < 0.05), ]
} else if (filt_by == "logFC"){
  sigtab <- res[abs(res$log2FoldChange) > 1, ]
} else if(filt_by == "all"){
  sigtab <- res
} else {
  stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
}
sigtab <- cbind(as(sigtab, "data.frame"),
                as(tax_table(phy)[rownames(sigtab), ], "matrix"))
# Log2FC plot
plot_dea(sigtab)
# Circos plot
## Prepare data
if (level == "Genus") {
  circus <- sigtab %>%
    select(Phylum, Genus, log2FoldChange) %>%
    dplyr::rename(up_level = Phylum, down_level = Genus) %>%
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
make_circos_plot(df = circus, group = links , col.pal = my_palette)

# Significance plot
df_plot <- phy %>%
  prune_taxa(taxa = rownames(sigtab)) %>%
  psmelt()
p = ggplot(df_plot, aes(x = eval(parse(text = level)), y = Abundance)) + 
  geom_boxplot(aes( fill = eval(parse(text = target))))+theme_light()+
  theme( axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank())+
  labs(fill = target)+
  facet_wrap( ~ eval(parse(text = level)), scales = 'free') +
  stat_compare_means(aes(group = eval(parse(text = target))), label = 'p.signif',vjust = 0.7, color = "red")
p = change_palette(p, viridis(n =length(unique(names(table(df_plot[target]))))))
p

sort(sigtab$Genus)

# require(phyloseq)
# require(tidyverse)
# library(broom)
# library(ggtext)
# phy = readRDS(file = "~/git/IBEROBDIA/02_preprocess/data/phy_Species.rds")
# df <- phy %>%
#   psmelt()
# 
# sig <- df %>%
#   nest(data = -Species) %>%
#   mutate(test = map(.x = data, ~ wilcox.test(Abundance ~ Statusv3, data = .x) %>% tidy)) %>%
#   unnest(test) %>%
#   filter(p.value < 0.05)
# 
# df %>%
#   inner_join(sig, by = 'Species') %>%
#   mutate(Species = str_replace(Species, "(.*)", "*\\1*"),
#          Statusv3 = factor(Statusv3, levels = c("No_DT2","DT2"))) %>%
#   ggplot(aes(x = Abundance, y = Species, color = Statusv3, fill = Statusv3)) +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.3,
#                                               dodge.width = 0.8 )) +
#   stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
#                geom = "pointrange",
#                position = position_dodge(width  = 0.8),
#                color = Statusv3, show.legend = FALSE) +
#   scale_x_log10() +
#   scale_color_manual(NULL,breaks = c('DT2', 'No_DT2'),
#                      values = viridis(3)[2:3],
#                      labels = c('DT2', 'No_DT2'))+
#   scale_fill_manual(NULL, values = viridis(3)[2:3],
#                     breaks = c('DT2', 'No_DT2'),
#                     labels= c('DT2', 'No_DT2')) +
#   labs(x = 'Abundance', y = NULL ) +
#   theme_classic() +
#   theme(
#     axis.text.y = element_markdown()
#   )



