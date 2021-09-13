#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --output=pvalue_check.out
#SBATCH --job-name=pvalue_check

Rscript pvalue_check.R