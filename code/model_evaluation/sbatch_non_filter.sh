#!/bin/bash
#SBATCH -J non_filter
#SBATCH --mem=10G
#SBATCH -t 24:00:00

module load miniconda
conda activate myPython
python /home/zc354/GRN/code/model_evaluation/non_filter_DF.py