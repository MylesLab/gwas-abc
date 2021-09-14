#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=ld_500_bins.out
#SBATCH --job-name=ld_500_bins

Rscript ld_500_bins.R