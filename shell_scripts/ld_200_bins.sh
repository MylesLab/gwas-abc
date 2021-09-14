#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=ld_200_bins.out
#SBATCH --job-name=ld_200_bins

Rscript ld_200_bins.R