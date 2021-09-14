#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --output=hetero_check.out
#SBATCH --job-name=hetero_check

Rscript hetero_check.R