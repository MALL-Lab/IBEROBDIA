# prepare data
require(microbiome)
extract_dataframe <- function(pseq){
  otu <- as.data.frame(t(as.data.frame(pseq@otu_table)))
  target <- pseq@sam_data$Status
  df <- cbind(otu, target)
  return(df)
}
level <- "Genus"
id <- "HvsPreDT2"
pseq <- readRDS(file = paste0("02_preprocess/data/pruned_", level,"_",id,"_pseq.rds"))
pseq <- subset_samples(physeq = pseq, Status != "DT2")
pseq <- microbiome::transform(x = pseq,transform = "clr")


df <- extract_dataframe(pseq = pseq)
write.csv2(x = df, file = paste0("03_SVM/data/df_",level,"_",
                                    id,".csv"),
           sep = ";", col.names = TRUE,row.names = TRUE)
table(df$target)
