setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lvl = "Species" #Species, Genus
if (lvl == "Genus") {
  exp = "first-experiment"
  out = "fourth-experiment"
}else if (lvl == "Species") {
  exp = "second-experiment"
  out = "fifth-experiment"
}else{
  stop("Wrong lvl introduced")
}
#paths
csv_path = paste0("~/git/IBEROBDIA/02_preprocess/data/df2csv/wout_outliers/")
rds_path = paste0("~/git/IBEROBDIA/03_training/toRun/",exp,"/")
out_path = paste0("~/git/IBEROBDIA/03_training/toRun/",out,"/")
#files
f.csv = list.files(path = csv_path, pattern = lvl)
f.rds = list.files(path = rds_path)
#lists
l.csv = list()
l.rds = list()

for (i in seq_along(f.csv)) {
  l.rds[[i]] = readRDS(file = paste(rds_path, f.rds[i], sep = "/"))
  l.csv[[i]] = read.csv2(file = paste0(csv_path, f.csv[i]), header = TRUE,sep = ";",row.names = 1)
  l.rds[[i]] = l.rds[[i]][rownames(l.csv[[i]]),]
  saveRDS(object = l.rds[[i]], file = paste0(out_path,f.rds[[i]]))
}

