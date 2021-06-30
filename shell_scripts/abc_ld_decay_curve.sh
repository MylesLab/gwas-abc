#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=abc_ld_decay_curve.out
#SBATCH --job-name=abc_ld_decay_curve

module load r
export R_LIBS=/project/6003429/myles_lab/bin/r_library

Rscript ld_decay_curve.R