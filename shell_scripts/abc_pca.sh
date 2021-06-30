#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=abc_pca.out
#SBATCH --job-name=abc_pca

module load java

/project/6003429/myles_lab/bin/tassel5/tasseladmin-tassel-5-standalone-68140c7ab0dc/run_pipeline.pl -Xmx40G -fork1 -plink -ped abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned.ped -map abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned.map -PrincipalComponentsPlugin -covariance true -ncomponents 10 -endPlugin -export abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned -runfork1