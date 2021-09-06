#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --output=background_ld.out
#SBATCH --job-name=background_ld

Rscript background_ld.R