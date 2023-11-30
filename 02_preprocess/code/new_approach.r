# Script to make DEA analysis 
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
level <- "Genus"
experiment <- "DT2_P_H"
levels <- c( "PreDT2", "Healthy", "DT2")
model <- '~Status + GS + BMI+ MS'
id <- "PreDT2vsDT2"
a <- 0.15
# Load phyloseq
phy <- readRDS(file = paste0("02_preprocess/data/phy_",
                             level,"_v2.rds"))

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
  p.adj.method = "fdr",
  prev.filter = 0, 
  mean.abund.filter = 0)
linda_result <- linda.res$output$StatusMS

df2plot <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj')) %>%
  rownames_to_column("names")

# Create a volcano  plot
ggplot(df2plot, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 2, color = ifelse((df2plot$log2FoldChange >= -1 & df2plot$log2FoldChange <= 1) & df2plot$padj > 0.2, "gray",
                                      ifelse(df2plot$padj <= 0.05, "red",
                                             ifelse(df2plot$padj <= 0.15, "orange",
                                                    ifelse(df2plot$padj <= 0.2, "yellow2", "green"))))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.15), linetype = "dashed", color = "orange") +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed", color = "yellow2") +
  annotate("text", x = Inf, y = -log10(0.06), label = "P-value = 0.05", vjust = -0.5, hjust = 1, size = 5) +
  annotate("text", x = Inf, y = -log10(0.16), label = "P-value = 0.15", vjust = -0.5, hjust = 1, size = 5) +
  annotate("text", x = Inf, y = -log10(0.21), label = "P-value = 0.2", vjust = -0.5, hjust = 1, size = 5) +
  geom_text_repel(data = subset(df2plot, padj <= 0.15), aes(label = names), size = 4, vjust = 1, hjust = 1) +
  labs(x = "Log Fold Change", y = "-log10(p-value)",title = model)  +
  theme_minimal() + theme(axis.text=element_text(size=12),
                          axis.title=element_text(size=12))

# select criteri to exclude
sigtab <- linda_result %>%
  dplyr::select(c('log2FoldChange', 'pvalue', 'padj','reject')) %>%
  mutate(log2FoldChange = abs(log2FoldChange)) %>%
  filter(padj < 0.05) 

#sigtab <- linda_result %>%
# mutate(log2FoldChange = abs(log2FoldChange)) %>%
# filter(pvalue < 0.05) %>%
# filter(log2FoldChange > 1)
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
sigtab
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
# saveRDS(object = circus, file = "figures/data/MS_Species.rds") # Save name to make a palette
make_circos_plot(df = circus, group = links , col.pal = my_palette)


# Heatmap
# phy <-  microbiome::transform(phy, transform =  'clr',target = "OTU" shift = 1)
phy <- phy %>%
  prune_taxa(taxa = rownames(sigtab))
df <- phy@sam_data
ordered_df <- df[order(df$Status ,decreasing = FALSE), ]
phy <- microViz::ps_reorder(ps = phy, sample_order = rownames(ordered_df))

t_names <- data.frame(tax_table(phy))$Genus
taxa_names(phy) <- t_names

# Annotation text size
ans = 7
# Rows text size
rns = 7
# Title size
ts = 9
# Legend title size
lts = 8
# legend text size
wls = 7

mat1 <- t(as.matrix(as.data.frame(phy@otu_table)))
common_min = min(c(mat2))
common_max = max(c(mat2))
col_fun = circlize::colorRamp2(c(-3,
                                 (0),
                                 3),
                               c("cornflowerblue","white", "brown3"))
# prepare complex heatmap

ann <- data.frame(phy@sam_data$Status)
colnames(ann) <- c('Status')
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
                            annotation_legend_param = list(
                              'Status' =list(nrow = 1,
                                             legend_direction = "horizontal",
                                             title_gp = gpar(fontsize = lts), 
                                             labels_gp = gpar(fontsize = wls),
                                             grid_height = unit(2, "mm"), 
                                             grid_width = unit(2, "mm"))),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            gap = unit(0.5, 'mm'))

mat2 <- apply(mat1, 2, function(x) log2(x + 0.1))

h1 <- Heatmap(scale(t(mat2)), name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "Discovery",
              column_title_gp = grid::gpar(fontsize = ts),
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)
h1










a = readRDS("~/git/BV_Microbiome/extdata/PRJNA242473/c380/PRJNA242473_tax_table.rds")
View(a)

tt <- as_tibble(PRJNA294119_tax_table)
View(tt %>% dplyr::filter(!is.na(Species)))


m <- as.matrix(as.data.frame(phy@otu_table))
# ALDEX2 
# matrix (rfeatures x csamples)

# Generate Monte Carlo samples of the Dirichlet distribution for each sample.
# Convert each instance using the centered log-ratio transform.
# This is the input for all further analyses.
set.seed(254)
x <- aldex.clr(m, phy@sam_data$Status)
# calculates expected values of the Welch's t-test and Wilcoxon rank
# test on the data returned by aldex.clr
x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)
# Determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)

par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.15
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.15
)

rownames_to_column(aldex_out, "genus") %>%
  filter(wi.eBH <= 0.25)  %>% # here we chose the wilcoxon output rather than tt
  dplyr::select(genus, we.eBH, wi.eBH, effect, overlap) %>%
  kable()


# ANCOMBC
#phyloseq
out <- ancombc2(
  data = phy,
  tax_level = "Species", 
  fix_formula = "Status", 
  p_adj_method = "fdr", 
  lib_cut = 0,
  prv_cut = 0,
  group = "Status", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.15, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
# store the FDR adjusted results [test on v1.2.2]
ancombc_result <- out$res %>%
  dplyr::select(starts_with(c("taxon", "lfc", "q")))
kable(ancombc_result)


# Maaslin2
# matrix (rsamples x cfeatures)
maaslin2_out <- Maaslin2(
  t(m),
  data.frame(phy@sam_data),
  output = "DAA example",
  transform = "AST",
  fixed_effects = model,
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "Status,MS",  
  normalization = "TSS",
  correction = "BH",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)
kable(head(filter(maaslin2_out$results, qval <= 0.3)))


# Linda
# dataframe (rfeatures x csamples)
meta <- as.tibble(phy@sam_data)
linda.res <- linda(
  as.data.frame(m), 
  meta, 
  formula = '~ Status + IR + GS + BMI', 
  alpha = 0.15,
  p.adj.method = "fdr",
  prev.filter = 0, 
  mean.abund.filter = 0)
linda_out <- linda.res$output$StatusNo_MS
View(linda_out)
# to scan the table for genera where H0 could be rejected:
kable(filter(as.data.frame(linda_out),
             abs(log2FoldChange) > 1))


# ZicoSeq
# matrix (rfeatures x csamples)
set.seed(123)
meta <- as.data.frame(as.matrix(phy@sam_data))
zicoseq.obj <- GUniFrac::ZicoSeq(meta.dat = meta, 
                                 feature.dat = m,
                                 grp.name = 'Status',
                                 adj.name = "GS", 
                                 feature.dat.type = 'count',
                                 prev.filter = 0,
                                 perm.no = 999,
                                 mean.abund.filter = 0,
                                 max.abund.filter = 0,
                                 return.feature.dat = T)
zicoseq_out <- cbind.data.frame(p.raw=zicoseq.obj$p.raw,
                                p.adj.fdr=zicoseq.obj$p.adj.fdr) 
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))
## x-axis is the effect size: R2 * direction of coefficient
ZicoSeq.plot(ZicoSeq.obj = zicoseq.obj,
             meta.dat = meta,
             pvalue.type ='p.raw')


df2plot <- linda_result %>%
  select(c( 'log2FoldChange', 'pvalue', 'padj'))
lfccut <- 1
p1 = EnhancedVolcano(toptable = df2plot,
                     lab = rownames(df2plot),
                     x ='log2FoldChange',
                     y = 'padj',
                     pCutoff = a,
                     FCcutoff = lfccut,
                     pointSize = 3.0,
                     labSize = 3.0,
                     ylim = c(0, -log10(10e-4)),
                     col = c("grey30", "forestgreen", "royalblue", "red2"),
                     colAlpha = 0.9,
                     cutoffLineType = 'twodash',
                     legendPosition = 'bottom',
                     title = paste0(experiment,"-", level, " vulcano"),
                     subtitle = model,
                     caption = paste("Log fold change cutoff,", lfccut,"; p-value cutoff,", a)
) 
p1