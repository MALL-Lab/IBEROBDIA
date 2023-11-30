# Script to make DEA analysis 
# -------------------------------------
# 3 variables:
# level : "Genus" or "Species"
# experiment: IMC SM Estado_Glucosa Glucosa_IMC SM_IMC G_SM_IMC Insulina_IMC
# filt_by: "pvalue" "logFC" "padj"
# -------------------------------------
# Declare variables 
level <- "Genus"
target <- "SM"
filt_by <- "padj"

# Required Packages 
require(phyloseq)
require(DESeq2)
require(ggpubr)
require(tidyverse)

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

# Load phyloseq
pseq <- readRDS(file = paste0("~/git/IBEROBDIA/02_preprocess/data/phy_",
                             level,".rds"))
df <-  read.csv2(file = "~/git/IBEROBDIA/new_data.csv", sep = ";", row.names = 1,
               header = TRUE)
m <-  data.frame(sample_data(pseq))
# Variables reorganization
samples = setdiff(x = rownames(df), y = rownames(m))
df = df[!(row.names(df) %in% samples),]
colnames(df)[2] <- "Estado_Glucosa"
colnames(df)[3] <- "IMC"
df = df[match(rownames(m), rownames(df)), ]
# Variables recategorization
df$SM[df$SM == "Sí"]<- 'SM'
df$SM[df$SM == "No"]<- 'No_SM'
df$Estado_Glucosa[df$Estado_Glucosa == "Sí"]<- 'GA'
df$Estado_Glucosa[df$Estado_Glucosa == "No"]<- 'GN'
df$IMC[df$IMC == "Sobrepeso 2"]<- 'OB'
df$IMC[df$IMC == "Obesidad "]<- 'OB'
df$IMC[df$IMC == "Normopeso"]<- 'NP'
df$Glucosa_IMC <- paste0(df$IMC,"_", df$Estado_Glucosa)
df$SM_IMC <- paste0(df$IMC,"_", df$SM)
df$G_SM_IMC <- paste0(df$IMC, "_", df$Estado_Glucosa, "_", df$SM)
colnames(df)[which(names(df) == "X")] <- "Insulina"
df$Insulina[df$Insulina == "Rango saludable" | df$Insulina == "Rango Saludable"] <- 'Sano'
df$Insulina[df$Insulina == "Resistencia significativa"] <- 'RS'
df$Insulina[df$Insulina == "Resistencia temprana" | df$Insulina == "Resistencia Temprana"] <- 'RT'
df$Insulina_IMC <- paste0(df$IMC,"_", df$Insulina)
m = cbind(m, df)
sample_data(pseq) <- m

if (target == "IMC") {
  pseq@sam_data$IMC = as.factor(pseq@sam_data$IMC)
  diaggen <- phyloseq_to_deseq2(pseq,  ~IMC)
  diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
  diaggen <- DESeq(diaggen, test = "Wald", fitType = "parametric")
  res <- results(diaggen, pAdjustMethod = "BH", cooksCutoff = FALSE)
  if (filt_by == "pvalue"){
    sigtab = res[which(res$pvalue < 0.05), ]
  } else if (filt_by == "padj"){
    sigtab = res[which(res$padj < 0.05), ]
  } else if (filt_by == "logFC"){
    sigtab <- res[abs(res$log2FoldChange) > 1, ]
  } else{
    stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
  }
  sigtab <- cbind(as(sigtab, "data.frame"),
                 as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
  taxas <- rownames(sigtab)
  # taxas <- c("ASV77", "ASV487","ASV541", "ASV525", "ASV223", "ASV763")
} else if (target == "SM"){
  pseq@sam_data$SM = as.factor(pseq@sam_data$SM)
  diaggen <- phyloseq_to_deseq2(pseq,  ~SM)
  diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
  diaggen <- DESeq(diaggen, test = "Wald", fitType = "parametric")
  res <- results(diaggen, pAdjustMethod = "BH", cooksCutoff = FALSE)
  if (filt_by == "pvalue"){
    sigtab = res[which(res$pvalue < 0.05), ]
  } else if (filt_by == "padj"){
    sigtab = res[which(res$padj < 0.05), ]
  } else if (filt_by == "logFC"){
    sigtab <- res[abs(res$log2FoldChange) > 1, ]
  } else{
    stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
  }
  sigtab <- cbind(as(sigtab, "data.frame"),
                  as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
  taxas <- rownames(sigtab)
  # taxas <- c("ASV6", "ASV163","ASV132", "ASV24", "ASV525", "ASV183","ASV162", "ASV2372", "ASV17")
} else if (target == "Estado_Glucosa"){
  pseq@sam_data$Estado_Glucosa = as.factor(pseq@sam_data$Estado_Glucosa)
  diaggen <- phyloseq_to_deseq2(pseq,  ~Estado_Glucosa)
  diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
  diaggen <- DESeq(diaggen, test = "Wald", fitType = "parametric")
  res <- results(diaggen, pAdjustMethod = "BH", cooksCutoff = FALSE)
  if (filt_by == "pvalue"){
    sigtab = res[which(res$pvalue < 0.05), ]
  } else if (filt_by == "padj"){
    sigtab = res[which(res$padj < 0.05), ]
  } else if (filt_by == "logFC"){
    sigtab <- res[abs(res$log2FoldChange) > 1, ]
  } else{
    stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
  }
  sigtab <- cbind(as(sigtab, "data.frame"),
                  as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
  taxas <- rownames(sigtab)
  # taxas <- c("ASV1895", "ASV132","ASV525", "ASV1938")
  
} else if (target == "Glucosa_IMC"){
  pseq@sam_data$Glucosa_IMC = as.factor(pseq@sam_data$Glucosa_IMC)
  diaggen <- phyloseq_to_deseq2(pseq,  ~Glucosa_IMC)
  diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
  diaggen <- DESeq(diaggen, test = "Wald", fitType = "parametric")
  res <- results(diaggen, pAdjustMethod = "BH", cooksCutoff = FALSE)
  if (filt_by == "pvalue"){
    sigtab = res[which(res$pvalue < 0.05), ]
  } else if (filt_by == "padj"){
    sigtab = res[which(res$padj < 0.05), ]
  } else if (filt_by == "logFC"){
    sigtab <- res[abs(res$log2FoldChange) > 1, ]
  } else{
    stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
  }
  sigtab <- cbind(as(sigtab, "data.frame"),
                  as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
  taxas <- rownames(sigtab)
  # taxas <- c("ASV487", "ASV1895","ASV132",
} else if (target == "SM_IMC") {
  pseq <- pseq %>%
    subset_samples(G_SM_IMC != "NP_GN_SM") 
  pseq@sam_data$SM_IMC = as.factor(pseq@sam_data$SM_IMC)
  diaggen <- phyloseq_to_deseq2(pseq,  ~SM_IMC)
  diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
  diaggen <- DESeq(diaggen, test = "Wald", fitType = "parametric")
  res <- results(diaggen, pAdjustMethod = "BH", cooksCutoff = FALSE)
  if (filt_by == "pvalue"){
    sigtab = res[which(res$pvalue < 0.05), ]
  } else if (filt_by == "padj"){
    sigtab = res[which(res$padj < 0.05), ]
  } else if (filt_by == "logFC"){
    sigtab <- res[abs(res$log2FoldChange) > 1, ]
  } else{
    stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
  }
  sigtab <- cbind(as(sigtab, "data.frame"),
                  as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
  taxas <- rownames(sigtab)
  # taxas <- c("ASV163", "ASV487", "ASV2850", "ASV6", "ASV24", "ASV77","ASV541",
  #          "ASV525","ASV2372", "ASV55")
  
} else if (target == "G_SM_IMC"){
  pseq <- pseq %>%
    subset_samples(G_SM_IMC != "NP_GN_SM")
  pseq@sam_data$G_SM_IMC = as.factor(pseq@sam_data$G_SM_IMC)
  diaggen <- phyloseq_to_deseq2(pseq,  ~G_SM_IMC)
  diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
  diaggen <- DESeq(diaggen, test = "Wald", fitType = "parametric")
  res <- results(diaggen, pAdjustMethod = "BH", cooksCutoff = FALSE)
  if (filt_by == "pvalue"){
    sigtab = res[which(res$pvalue < 0.05), ]
  } else if (filt_by == "padj"){
    sigtab = res[which(res$padj < 0.05), ]
  } else if (filt_by == "logFC"){
    sigtab <- res[abs(res$log2FoldChange) > 1, ]
  } else{
    stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
  }
  sigtab <- cbind(as(sigtab, "data.frame"),
                  as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
  taxas <- rownames(sigtab)
  # taxas <- c("ASV163", "ASV6","ASV525", "ASV55")
} else if (target == "Insulina_IMC"){
  pseq <- pseq %>%
    subset_samples(Insulina != "Null") %>%
    subset_samples(Insulina_IMC != "NP_RS") %>%
    subset_samples(Insulina_IMC != "NP_RT")
  pseq@sam_data$Insulina_IMC = as.factor(pseq@sam_data$Insulina_IMC)
  diaggen <- phyloseq_to_deseq2(pseq,  ~Insulina_IMC)
  diaggen <- estimateSizeFactors(diaggen, type = 'poscounts')
  diaggen <- DESeq(diaggen, test = "Wald", fitType = "parametric")
  res <- results(diaggen, pAdjustMethod = "BH", cooksCutoff = FALSE)
  if (filt_by == "pvalue"){
    sigtab = res[which(res$pvalue < 0.05), ]
  } else if (filt_by == "padj"){
    sigtab = res[which(res$padj < 0.05), ]
  } else if (filt_by == "logFC"){
    sigtab <- res[abs(res$log2FoldChange) > 1, ]
  } else{
    stop("Please insert a valid filt parameter (pvalue, padj or logFC)")
  }
  sigtab <- cbind(as(sigtab, "data.frame"),
                  as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
  taxas <- rownames(sigtab)
  #taxas = c("ASV413", "ASV223")
} else {
  stop("Please insert a valid target")
}

# Log2FC plot
plot_dea(sigtab)
# Circos plot
## Prepare data
if (level == "Genus") {
  circus <- sigtab %>%
    select(Phylum, Genus, log2FoldChange) %>%
    rename(up_level = Phylum, down_level = Genus) %>%
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
    rename(up_level = Genus, down_level = Species) %>%
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
df_plot <- pseq %>%
  prune_taxa(taxa = taxas) %>%
  psmelt()
p <- ggplot(df_plot, aes(x = eval(parse(text = level)), y = Abundance)) + 
  geom_boxplot(aes( fill = eval(parse(text = target))))+theme_light()+
  theme( axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank())+
  labs(fill = target)+
  facet_wrap( ~ eval(parse(text = level)), scales = 'free') +
  stat_compare_means(aes(group = eval(parse(text = target))), label = 'p.signif',vjust = 0.7)
p <- change_palette(p, viridis(n =length(unique(names(table(df_plot[target]))))))
p


