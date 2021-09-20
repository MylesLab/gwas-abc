#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --output=alt_tassel_manhattan.out
#SBATCH --job-name=alt_tassel_manhattan
#SBATCH --mail-user=sophie.watts@dal.ca

Rscript alt_tassel_manhattan.R