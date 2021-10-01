#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --output=phenolics_mlmm.out
#SBATCH --job-name=phenolics_mlmm

Rscript phenolics_mlmm.R