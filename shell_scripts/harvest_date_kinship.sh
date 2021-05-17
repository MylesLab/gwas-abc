#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --output=harvest_date_kinship.out
#SBATCH --job-name=harvest_date_kinship

module load java

/project/6003429/myles_lab/bin/tassel5/tasseladmin-tassel-5-standalone-68140c7ab0dc/run_pipeline.pl -plink -ped abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001.ped -map abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001.map -KinshipPlugin -method Centered_IBS -endPlugin -export abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001_kinship.txt -exportType SqrMatrix