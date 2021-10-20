#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=harvest_date_pca.out
#SBATCH --job-name=harvest_date_pca

module load java

/project/6003429/myles_lab/bin/tassel5/tasseladmin-tassel-5-standalone-68140c7ab0dc/run_pipeline.pl -Xms512m -Xmx10g -fork1 -plink -ped abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001.ped -map abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001.map -PrincipalComponentsPlugin -covariance true -endPlugin -export abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001_ -runfork1