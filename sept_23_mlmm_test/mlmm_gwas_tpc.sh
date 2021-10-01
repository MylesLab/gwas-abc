#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --mem=64G
#SBATCH --output=mlmm_gwas_tpc.out
#SBATCH --job-name=mlmm_gwas_tpc

Rscript mlmm_gwas_tpc.R