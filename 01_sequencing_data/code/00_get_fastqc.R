#FASTQC
# 16S Fastq Processing
library(dada2)
path = "projects/IBEROBDIA/ARCHIVOS_FASTQ"
fnFs <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))

## Gonna check fasta quality with Fastqc 
l = list()
# make command line for each raw
for (i in seq_along(fnFs)){
  l[[i]] = paste("Fastqc", fnFs[i])
}
# Run those commands commmand 
for (i in seq_along(fnFs)){
  system(l[[i]])
}

# Prepare open htmls informs from Fastqc
htmls <- sort(list.files(path, pattern=".html", full.names = TRUE))
open_htmls = htmls
for (i in seq_along(htmls)){
  open_htmls[i] = paste("open", open_htmls[i])
}
system(open_htmls[1])
for (i in seq.int(20,30)) {
  system(open_htmls[i])
}
