#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=64G
#SBATCH --output=mlmm_gwas_acidity.out
#SBATCH --job-name=mlmm_gwas_acidity

Rscript mlmm_gwas_acidity.R