#!/bin/bash

#SBATCH -t 0-96:00:00
#SBATCH --mem=80G
#SBATCH -p general,pi_zhao
#SBATCH -J EPV_TF_unnomalization

module load miniconda
conda activate myPython
python3 /gpfs/gibbs/pi/zhao/zc354/GRN/code/Early_precision_value_new.py