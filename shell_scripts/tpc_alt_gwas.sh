#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --output=tpc_alt_gwas.out
#SBATCH --job-name=tpc_alt_gwas

module load java

/project/6003429/myles_lab/bin/tassel5/tasseladmin-tassel-5-standalone-68140c7ab0dc/run_pipeline.pl -fork1 -plink -ped /project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/tpc_geno_filtered.ped -map /project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/tpc_geno_filtered.map -fork2 -r /project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/phenotype_data/tassel_format/tpc_pheno_reformated.txt —fork3 -importGuess /project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/pca/tpc_pca_1.txt -combine4 -input1 -input2 -input3 -intersect -fork5 -importGuess /project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/kinship/tpc_kinship.txt -combine6 -input4 -input5 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None —mlmOutputFile /project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/alt_tassel_gwas/tpc