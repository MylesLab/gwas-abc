#!/usr/bin/bash

#SBATCH --ntasks=21
#SBATCH --mem-per-cpu=5G
#SBATCH --account=def-smyles
#SBATCH --job-name=run-mlmm-rscripts
#SBATCH --output=run-mlmm-rscripts.out
#SBATCH --time=3-00:00:00

#SBATCH --mail-user=sophie.watts@dal.ca
#SBATCH --mail-type=ALL

srun Rscript run_mlmm_gwas.R juiciness_16_harv
srun Rscript run_mlmm_gwas.R flowering_jul_16_harv
srun Rscript run_mlmm_gwas.R precocity_16_harv
srun Rscript run_mlmm_gwas.R acidity_17_stor
srun Rscript run_mlmm_gwas.R brix_17_stor
srun Rscript run_mlmm_gwas.R weight_avg_17_stor
srun Rscript run_mlmm_gwas.R firmness_avg_17_stor
srun Rscript run_mlmm_gwas.R brix_acid_17_stor
srun Rscript run_mlmm_gwas.R brix_17_harv
srun Rscript run_mlmm_gwas.R acidity_17_harv
srun Rscript run_mlmm_gwas.R date_jul_17_harv
srun Rscript run_mlmm_gwas.R firmness_avg_17_harv
srun Rscript run_mlmm_gwas.R weight_avg_17_harv
srun Rscript run_mlmm_gwas.R brix_acid_17_harv
srun Rscript run_mlmm_gwas.R percent_acidity_17
srun Rscript run_mlmm_gwas.R percent_brix_17
srun Rscript run_mlmm_gwas.R percent_weight_avg_17
srun Rscript run_mlmm_gwas.R percent_firmness_avg_17
srun Rscript run_mlmm_gwas.R percent_brix_acid_17
srun Rscript run_mlmm_gwas.R tpc
srun Rscript run_mlmm_gwas.R time_ripen_2017
