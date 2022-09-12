# Load packages
library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)

# Load phy.object
phy = readRDS("projects/IBEROBDIA/Data/TL_251/TL_251_phy.rds")

# Summary
ntaxa(phy)
nsamples(phy)
sample_names(phy)[1:5]  
rank_names(phy)  
sample_variables(phy)  
otu_table(phy)[1:5, 1:5]  
tax_table(phy)[1:5, 1:4]

p = data.frame(sample_data(phy))
# R adding a column to dataframe based on values in other columns:
p<- p %>% 
  mutate(New_Status = if_else(Status == "Health", "Health", "Affected"))
sample_data(phy) <- p

# Preprocess Phyloseq

#1. Aglomerar en genero (Hay muy pocas especies, si aglomeramos por estas se va la mayoria al carajo)
g_phy = tax_glom(physeq = phy,taxrank=rank_names(phy)[7])

#2. Filtrar aquellas ASVs que no cumplan unos criterios
## Make and apply the filter
filter <- phyloseq::genefilter_sample(g_phy, filterfun_sample(function(x) x > 0), 
                                      A = 0.05*nsamples(g_phy))
g_phy <- prune_taxa(filter, g_phy)
## If we want filter by variance 
g_phy = filter_taxa(g_phy, function(x) var(x) > 0.1, TRUE)


#3. Sacar PLOTS alfa,beta,HM (Ordenar por target), Recatalogar y mirar la disposicion de generos/especies Health vs all 
library(microbiome)
library(ggplot2)
#Reads Dsitribution
sample_sum_df <- data.frame(sum = sample_sums(g_phy))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) 
#

df <- psmelt(g_phy)

#4. ML con NP vs OB, Health vs All??
plot_richness(phy, color = "Status", x = "CI", measures = c("Observed", "Chao1", "Shannon")) + geom_boxplot(aes(fill = Status), alpha=.7)

ord.nmds.bray <- ordinate(phy, method="NMDS", distance="bray")
plot_ordination(phy, ord.nmds.bray, color="CI", title="Bray NMDS")


# Necesitamos obtener las taxa más abundantes, en este caso el top 15
library(fantaxtic)
library(microViz)
k = ps_arrange(g_phy, Status)
top15 <- get_top_taxa(physeq_obj = k, n = 20, relative = T,
                      discard_other = F, other_label = "Other")
# Ya que no todas las taxa fueron clasificadas a nivel de especie, generamos etiquetas compuestas de distintos rangos taxonómicos para el gráfico
top15 <- name_taxa(top15, label = "", species = T, other_label = "Other")
# Finalmente graficamos
fantaxtic_bar(top15, color_by = "Species", label_by = "Species", facet_by =NULL, grid_by = NULL,order_alg = "as.is", other_color = "Grey") -> ptop15
ptop15



