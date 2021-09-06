#!/bin/bash
#SBATCH --time=3-2:00:00
#SBATCH --mem=64G
#SBATCH --output=mlmm_gwas.out
#SBATCH --job-name=mlmm_gwas

Rscript mlmm_gwas.R