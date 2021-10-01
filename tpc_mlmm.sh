#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --output=tpc_mlmm.out
#SBATCH --job-name=tpc_mlmm

#SBATCH --mail-user=sophie.watts@dal.ca
#SBATCH --mail-type=ALL

Rscript tpc_mlmm.R