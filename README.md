# abc-gwas
Genome-wide association study on fruit quality and phenology traits from the ABC.

## Source
**abc_ploidy** script for looking at heterozygosity per individual to examine ploidy in the ABC.

**snps_stats** script for analyzing SNP stats.

**ld_decay** R script run to caluclate LD decay.

**pheno_prep** R script to create phenotype files and apple id files for running GWAS.

**pop_analyses** script for visualizing PCA and LD decay.

**phenotype_curation** script for prepping phenotype table.

## Shell Scripts

**abc_ld_decay_curve** script to run ld_decay.R

**abc_pca** script to run PCA with TASSEL.

## Outputs

**contig_snps_remove.txt** list of unanchored SNPs to be removed from SNP table.

**geno_pheno_meta_data.csv** Filtered phenotype table that includes apple IDs with phenotype data and genotype data.

**pheno_list.txt** list with names of phenotypes for GWAS.