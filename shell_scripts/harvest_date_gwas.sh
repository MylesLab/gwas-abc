#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=harvest_date_gwas.out
#SBATCH --job-name=harvest_date_gwas

module load r
export R_LIBS=/project/6003429/myles_lab/bin/r_library

Rscript harvest_date_gwas.R