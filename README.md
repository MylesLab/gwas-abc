# abc-gwas
Genome-wide association study on fruit quality and phenology traits from the ABC.

## Data

**41438\_2019\_190\_MOESM2\_ESM (2).xls** Phenolic data from Kendra's polyphenol manuscript.

**41438\_2019\_190\_MOESM3\_ESM (4).xls** Phenolic data from Kendra's polyphenol manuscript.

**vineland_snps** SNP names for 8 SNPs genotyped at Vineland.

**full\_pheno_names** Short and full names for the traits measured in the ABC.

**pheno_header** Text file with the header to be added to the phenoytpe data files to format them for TASSEL.

**pheno\_meta_data.csv** Phenotype data and meta data for all accessions in the ABC.

**abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001.frq** Minor allele frequency file for ABC snps.

**abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned1** Results from PCA on full ABC snp set (contig snps removed and LD pruned).

**abc\_combined\_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned2** Eigen values from PCA on full ABC snp set (contig snps removed and LD pruned).

**abc\_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001.nosex** File with the apple IDS that have genetic data.

**abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_het.het** Plink file with heterozygosity per individual.


## Outputs

**contig_snps_remove.txt** list of unanchored SNPs to be removed from SNP table.

**geno_pheno_meta_data.csv** Filtered phenotype table that includes apple IDs with phenotype data and genotype data.

**pheno_list.txt** list with names of phenotypes for GWAS.

**phenotypes4gwas** contains two files for each phenotype: one with a list of apple IDs for that phenotype that have trait data and the second a file with the trait measurements for that phenotype.

**phenotype_sample_sizes** table with samples sizes for each phenotype.

**simple_m.out** Output from Simple M package that calculates the effective number of markers.

**pc1\_vs\_phenos.csv** The R-square and p-values from the correlation of phenotypes with PC1.

**pc2\_vs\_phenos.csv** The R-square and p-values from the correlation of phenotypes with PC2.

## Shell Scripts

**ld_decay_curve.sh** script to run ld_decay.R

**background_ld.sh** script to back_ground_ld.R

**simple_m.sh** script to run simple_m.R

**abc_pca.sh** script to run PCA for the whole SNP set with TASSEL.

**genotype_filtering_plink.sh** script that contains PLINK commands to filter the ABC MAP and PED to only containing apple IDs for a particular phenotype and outputs a MAP and PED for each phenotype.

**kinship.sh** script that contains the tassel commands to make a kinship matrix for each phenotype.

**pca.sh** script commands to run PCA with tassel.

**tassel_gwas.sh** script with commands to run TASSEL GWAS for each phenotype.

**geno_raw.sh** script with commands to recode PED and MAP files into .raw files for the MLMM gwas.

**mlmm_gwas** script to run the MLMM GWAS code.

**phenolics_mlmm** script to run the code to test out mlmm for just phenolic content.

**pval_sorting** script to sort the p-values from TASSEL GWAS from most significant to least significant.


## Source
**abc_ploidy.Rmd** script for looking at heterozygosity per individual to examine ploidy in the ABC.

**snps_stats.R** script for analyzing SNP stats.

**ld_decay_curve.R** R script run to caluclate LD decay.

**pop_analyses.Rmd** script for visualizing PCA and LD decay.

**phenotype_curation.Rmd.Rmd** script for prepping phenotype table and outputting the list of phenotypes.

**background_ld.R** script to calculate background linkage disequilibrium decay.

**gwas_script.Rmd** script to create shell scripts to filter phenotype and  genotype data, create kinship matrices, run PCA, run TASSEL GWAS.

**mlmm_gwas.R** script to create shell scripts to run MLMM GWAS.

**pheno_reformat_tassel.R** script to reformat the phenotype files for TASSEL.

**simple_m.R** script to run Simple M to calculate the effective number of markers.

**snp_stats.R** script to analyze the SNP stats of ABC such as MAF distribution, number of SNPs etc.

**tassel_manhattans.R** script to write shell scripts to plot the output of the TASSEL GWAS as Manhattan plots and QQ plots.

**phenolics_deep_dive.R** script to compare phenolic content in the ABC to phenolic content measured by Kendra in the polyphenol paper.

**correlations.R** Script to run correlations between phenotypes.

**phenolics_mlmm.R** Code to test the mlmm package to run GWAS just for phenolics.

**pval_sorting.R** Code to sort p-values from TASSEL.





