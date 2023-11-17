#!/bin/bash

#SBATCH -t 0-72:00:00
#SBATCH --mem=40G
#SBATCH -p general,pi_zhao
#SBATCH -J gen_gene_peak_matrix

module load R/4.2.0-foss-2020b
Rscript ~/GRN/code/data_produce_linkpeaks.R