#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=simple_m.out
#SBATCH --job-name=simple_m

Rscript simple_m.R