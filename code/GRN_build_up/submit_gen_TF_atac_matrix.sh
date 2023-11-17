#!/bin/bash
#SBATCH -t 0-72:00:00
#SBATCH --mem=80G
#SBATCH -p general,pi_zhao
#SBATCH -J sct

module load R/4.2.0-foss-2020b
Rscript ~/GRN/code/bmmc_norm.R






