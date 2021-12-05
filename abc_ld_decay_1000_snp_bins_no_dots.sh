#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=abc_ld_decay_1000_snp_bins_no_dots.out
#SBATCH --job-name=abc_ld_decay_1000_snp_bins_no_dots

Rscript abc_ld_decay_1000_snp_bins_no_dots.R