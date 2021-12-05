#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=ld_decay_200.out
#SBATCH --job-name=ld_decay_200

Rscript ld_decay_200.R