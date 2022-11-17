

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

<div align="center">
    <h1>Genome-wide association study on fruit quality and phenology traits in Canada's Apple Biodiversity Collection</h1>
    <br />
    Sophie Watts<sup>1</sup>, ‪Zoë Migicovsky<sup>1</sup>, Gavin M. Douglas<sup>1</sup>, Sean Myles*<sup>1</sup>
    <br />
</div>
<b>Affiliations</b><br />
<sup>1</sup>Department of Plant, Food, and Environmental Sciences, Faculty of Agriculture, Dalhousie University

---

## Table of contents

1. [About](#1-about) </br>
2. [Abstract](#2-abstract) </br>
3. [Analysis](#4-analysis) </br>
3a. [File Architecture](#4a-file-architecture) <br />
3b. [Getting the code](#4b-getting-the-code) <br />
3c. [Dependencies](#4c-dependencies) <br />
3d. [Reproducing the results](#4d-reproducing-the-results) <br />

## 1. About
This repository contains data and scripts used to repoduce analyses in the manuscript "Genome-wide association studies in Canada's apple biodiversity collection".

## 2. Abstract

## 3. Analysis
### 3a. File architecture

```
├── data
├── source
├── outputs
└── figures
```

- The `data` directory contains all the raw data that was used for this project.
- The `outputs` directory contains the files generated from the raw data through data curation.
- The `source` directory contains the scripts used for data curation and performing various analyses of this project.
- The `figures` directory contains the intermediary figures generated for this project. Final figures were assembled using Adobe Illustrator. 


### 3b. Getting the code

You can download a copy of all the files in this repository by cloning the git repository:

```sh
$ git clone https://github.com/MylesLab/abc-gwas.git
```

### 3c. Dependencies

The library dependencies for running the R code are below:

```
tidyverse
viridis

```

### 3d. Reproducing the results

Given that you have above dependencies installed in your working machine, you should be able to run any of the code in this repository. However, in order to reproduce the analyses and obtain the same results, here is the list and order of the scripts you should run:

**1. Data Curation**
<br />
You first need to generate the dataset which is going to be used for running the analyses.<br />
**1a.** `data_curation/generate-final-dessert-apples-list.R` - Running this file will generate `data/processed/final_dessert_apple_phenotype_data.tsv`.<br />
**1b.** `data_curation/generate-final-cider-apples-list.R` - Running this file will generate `data/processed/final_cider_apple_phenotype_data.tsv`. <br />
**1c.** `data_curation/generate-final-phenotype-table.R` - Running this file will generate `data/processed/final_phenotype_table.tsv` which is the main file used for all the analyses. <br />

**2. Analyses** 
<br />
Once the `data/processed/final_phenotype_table.tsv` file is generated, we can start to run the anlyses. <br />
**2a.** `analyses/pca-density-analysis.R` - Running this file will generate Figure-2_density plot, as well as some of the intermediary data files which are required for generating final Figure 1 plot. <br />
**2b.** `analyses/distance-analysis.R` - Running this file will generate the final Figure 1 plot <br />




