setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(phyloseq)
library(ggfortify)
library(ggplot2)
library(viridis)
library(plotly)
lvl = "Species" # Genus or Species
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
ggplot2::autoplot(pca.re, data = df, colour = "Statusv3",
                  shape = "Statusv3",size = 3)


var = apply(df[5:ncol(df)][,3:length(colnames(df[5:ncol(df)]))],
            2, function(x) var(na.omit(x)))
var.ordenado <- sort(var,decreasing = T)
names(var.ordenado[1:3])

pca3D = plot_ly(x = df$Prevotella_9_copri, 
                   y = df$Bacteroides_uniformis, 
                   z = df$Bacteroides_caccae,
                   type = "scatter3d", 
                   mode = "markers",colors =  viridis(2),
                   color = as.factor(df$Statusv3)) %>%
  layout(legend = list(size = 35,orientation = "h", x = 0.3, y = 0),
         scene = list(xaxis = list(title = "Prevotella_9_copri", tickfont = list(size = 12)),
                      yaxis = list(title = "Bacteroides_uniformis", tickfont = list(size = 12)),
                      zaxis = list(title = "Bacteroides_caccae", tickfont = list(size = 12))))
pca3D


####
df = readRDS(paste0("~/git/IBEROBDIA/02_preprocess/data/df_",lvl,".rds"))
rmv = c("Status", "Gender", "Statusv2", "Statusv3")
df <- df[, ! names(df) %in% rmv, drop = F]
names(df)[names(df) == "CI"] <- "target"
saveRDS(object = df, paste0("~/git/IBEROBDIA/03_training/toRun/DF_CI_",lvl,".rds"))

df = readRDS(paste0("~/git/IBEROBDIA/02_preprocess/data/df_",lvl,".rds"))
rmv = c("CI", "Gender", "Statusv2", "Statusv3")
df <- df[, ! names(df) %in% rmv, drop = F]
names(df)[names(df) == "Status"] <- "target"
saveRDS(object = df, paste0("~/git/IBEROBDIA/03_training/toRun/DF_Status_",lvl,".rds"))

df = readRDS(paste0("~/git/IBEROBDIA/02_preprocess/data/df_",lvl,".rds"))
rmv = c("CI", "Gender", "Status", "Statusv3")
df <- df[, ! names(df) %in% rmv, drop = F]
names(df)[names(df) == "Statusv2"] <- "target"
saveRDS(object = df, paste0("~/git/IBEROBDIA/03_training/toRun/DF_Healthy_",lvl,".rds"))

df = readRDS(paste0("~/git/IBEROBDIA/02_preprocess/data/df_",lvl,".rds"))
rmv = c("CI", "Gender", "Status", "Statusv2")
df <- df[, ! names(df) %in% rmv, drop = F]
names(df)[names(df) == "Statusv3"] <- "target"
saveRDS(object = df, paste0("~/git/IBEROBDIA/03_training/toRun/DF_DT2_",lvl,".rds"))


####
####
####
####
