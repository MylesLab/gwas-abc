#!/usr/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --account=def-smyles
#SBATCH --time=5:00:00

Rscript mlmm_gwas_batch_script_final.R $1

