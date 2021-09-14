#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=ld_100_bins.out
#SBATCH --job-name=ld_100_bins

Rscript ld_100_bins.R