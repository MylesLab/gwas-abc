#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=maf_pheno_cor.out
#SBATCH --job-name=maf_pheno_cor

Rscript maf_pheno_cor.R