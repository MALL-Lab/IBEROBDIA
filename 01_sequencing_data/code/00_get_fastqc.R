# FASTQC
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
## 16S Fastq quality check
path = "../../extdata/FASTQ"
fnFs <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))

## Gonna check fasta quality with Fastqc 
l = list()
# make command line for each raw (need FASTQC installed and ad to envir)
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
for (i in seq.int(1,2)) {
  system(open_htmls[i])
}
