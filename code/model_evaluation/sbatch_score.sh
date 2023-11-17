#!/bin/bash
#SBATCH -J score
#SBATCH --mem=10G
#SBATCH -t 48:00:00

module load miniconda
conda activate myPython
python /home/zc354/GRN/code/model_evaluation/score_DF.py