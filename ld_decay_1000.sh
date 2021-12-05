#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=ld_decay_1000.out
#SBATCH --job-name=ld_decay_1000

Rscript ld_decay_1000.R