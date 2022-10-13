#!/bin/bash
#SBATCH -N 1
#SBATCH -p thinnodes
#SBATCH -t 01:00:00
#SBATCH --mem=64G
#SBATCH --error="/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/IBEROBDIA/01_sequencing_data/out/16S_error.txt"
#SBATCH --output="/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/IBEROBDIA/01_sequencing_data/out/16S_salida.txt"

module load gcc/6.4.0 R/4.0.2
TMP=/mnt/netapp2/Store_uni/home/ulc/co/dfe/TEMPR Rscript /mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/01_sequencing_data/code/01_pipeline_16S.r