#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=ld_decay_curve.out
#SBATCH --job-name=ld_decay_curve

Rscript ld_decay_curve.R