setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(phyloseq)
library(ggfortify)
library(ggplot2)
library(viridis)
library(plotly)
lvl = "Genus" # Genus or Species
phy = readRDS(paste0("../data/phy_",lvl ,".rds"))

#1. PCA
#First of all we extract the data fram object with genus or species names instead of ASVs
otu = as.data.frame(t(phy@otu_table@.Data))
tax = as.data.frame(phy@tax_table@.Data)
tax$otu_names <- rownames(tax)
tax = as.data.frame(cbind(otus = tax$otu_names, genus = tax[,lvl]))
#identical(rownames(otu), tax$otus)
rownames(otu) = tax$genus
meta = as.data.frame(phy@sam_data)
#identical(rownames(meta), colnames(otu))
otu = as.data.frame(t(otu))
df = as.data.frame(cbind(meta, otu))
# Saving that dataframe
saveRDS(object = df, file = paste0("../data/df_",lvl,".rds"))

pca.re <- stats::prcomp(df[6:ncol(df)],
                        center = TRUE,
                        scale. = FALSE)
ggplot2::autoplot(pca.re, data = df, colour = "Status",
                  shape = "Status",size = 3)


var = apply(df[5:ncol(df)][,3:length(colnames(df[5:ncol(df)]))],
            2, function(x) var(na.omit(x)))
var.ordenado <- sort(var,decreasing = T)
names(var.ordenado[1:3])

pca3D = plot_ly(x = df$Prevotella_9, 
                   y = df$Bifidobacterium, 
                   z = df$`UCG-002`,
                   type = "scatter3d", 
                   mode = "markers",colors =  viridis(4),
                   color = as.factor(df$Status)) %>%
  layout(legend = list(size = 35,orientation = "h", x = 0.3, y = 0),
         scene = list(xaxis = list(title = "Prevotella_9", tickfont = list(size = 12)),
                      yaxis = list(title = "Bifidobacterium", tickfont = list(size = 12)),
                      zaxis = list(title = "UCG-002", tickfont = list(size = 12))))
pca3D


####
df = readRDS("git/IBEROBDIA/02_preprocess/data/df_Genus.rds")
rmv = c("Status", "Gender", "Statusv2", "Statusv3")
df <- df[, ! names(df) %in% rmv, drop = F]
names(df)[names(df) == "CI"] <- "target"
df

saveRDS(object = df, "git/IBEROBDIA/03_training/toRun/DF_CI_genus.rds")

####
####
####
####
