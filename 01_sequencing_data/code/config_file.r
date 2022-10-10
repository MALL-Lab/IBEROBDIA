# DADA2 config file
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/IBEROBDIA')

experiment.name = "TL_251"
input.dir.path =  "extdata/FASTQ"
out.dir.path  = "01_sequencing_data/data"

out.path = paste(out.dir.path, experiment.name, '/', sep = '')

if (dir.exists(out.path) == FALSE) {
  dir.create(out.path)
  message(paste('Creating ', ExperimentName, ' directory!'))
}

# Filter parameters
dada2.truncLen = 251
dada2.maxN = 0
dada2.maxEE = 3
dada2.truncQ = 2
dada2.trimLeft = 15

# Learn Errors parameters
dada2.nbases = 1e+9

# DADA2 core algorithm parameters
dada2.homopolymer_gap_penalty = -1
dada2.band_size = 32

# Removing Quimeras parameters
dada2.method = "consensus"

# Taxonomical assignment
dada2.path.refFasta = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/tax/silva_nr99_v138.1_train_set.fa.gz"
dada2.tryRC = TRUE
dada2.path.refFasta.species = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/tax/silva_species_assignment_v138.1.fa.gz"
dada2.tryRC.species = FALSE


