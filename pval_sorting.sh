#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --output=pval_sorting.out
#SBATCH --job-name=pval_sorting

Rscript pval_sorting.R