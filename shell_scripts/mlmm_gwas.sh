#!/bin/bash
#SBATCH --time=3-6:00:00
#SBATCH --mem=64G
#SBATCH --output=mlmm_gwas.out
#SBATCH --job-name=mlmm_gwas

Rscript mlmm_gwas.R