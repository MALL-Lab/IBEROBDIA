# FASTQC
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
## 16S Fastq quality check
path = "../../extdata/FASTQ"
fnFs <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))

## Gonna check fasta quality with Fastqc 
l = list()
# make command line for each raw (need FASTQC installed and adsd alias to env)
for (i in seq_along(fnFs)){
  l[[i]] = paste("Fastqc", fnFs[i])
}
# Run those commands commmand 
for (i in seq_along(fnFs)){
  system(l[[i]])
}
