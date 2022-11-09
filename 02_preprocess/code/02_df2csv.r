setwd(dirname(rstudioapi::getSourceEditorContext()$path))
exp = "second-experiment" #first experiment, second-experiment
pth = paste0("~/git/IBEROBDIA/03_training/toRun/",exp)
l = list.files(path = pth)
l.df = list()
for (i in seq_along(l)) {
l.df[[i]] = readRDS(paste(pth,l[i],sep = "/"))
write.csv2(x =l.df[[i]], file = paste0("../data/df2csv/",substr(x =l[i], start = 1, stop =  nchar(l[i])-4), ".csv"))
}
