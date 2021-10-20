#!/usr/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --account=def-smyles
#SBATCH --time=2:00:00

#SBATCH --mail-user=sophie.watts@dal.ca
#SBATCH --mail-type=ALL

Rscript mlmm_gwas_batch_script_jpeg.R $1
