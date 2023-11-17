#!/bin/bash

#SBATCH -t 0-24:00:00
#SBATCH --mem=40G
#SBATCH -p general,pi_zhao
#SBATCH -J distance

module load R/4.2.0-foss-2020b
Rscript /gpfs/gibbs/pi/zhao/zc354/GRN/code/trio/cal_500kb_all.R