<div align="center">
    <h1>Genome-wide association studies in Canada's Apple Biodiversity Collection</h1>
    <br />
    Sophie Watts<sup>1</sup>, ‪Zoë Migicovsky<sup>1</sup>, Sean Myles*<sup>1</sup>
    <br />
</div>
<b>Affiliations</b><br />
<sup>1</sup>Department of Plant, Food, and Environmental Sciences, Faculty of Agriculture, Dalhousie University

---

## Table of contents

1. [About](#1-about) </br>
2. [Abstract](#2-abstract) </br>
3. [Analysis](#3-analysis) </br>
3a. [Getting the code](#3a-getting-the-code) <br />
3b. [File architecture](#3b-file-architecture) <br />
3c. [Description of files](#3c-description-of-files) <br />

## 1. About
This repository contains data and scripts used to repoduce analyses in the manuscript "Genome-wide association studies in Canada's apple biodiversity collection".

## 2. Abstract
Apple fruit quality traits such as fruit texture, sugar content, and firmness retention during storage are key targets for breeders. Understanding the genetic control of fruit quality traits can enable the development of genetic markers, useful for marker-assisted breeding of new apple cultivars. We genotyped over 260,000 single nucleotide polymorphisms (SNPs) across 1,054 apple accessions from Canada’s Apple Biodiversity Collection and performed genome-wide association for 21 fruit quality and phenology traits. We identified a locus on chromosome 15 associated with phenolic content and a locus on chromosome 10 that is associated with softening. We demonstrate that the top SNP on chromosome 10 is a better predictor of softening than markers commonly used for marker-assisted breeding of this trait. In addition, we identified a single locus on chromosome 3 that is associated with numerous traits including ripening time, firmness at harvest, and firmness after storage. The top SNP at the chromosome 3 locus is a nonsynonymous mutation within the NAC18.1 transcription factor. Given the association between variation at NAC18.1 and several key traits, we propose a model for the allelic effects at NAC18.1 on apple ripening and softening.


## 3. Analysis

### 3a. Getting the code

You can download a copy of all the files in this repository by cloning the git repository:

```sh
$ git clone https://github.com/MylesLab/abc-gwas.git
```


### 3b. File architecture

```
├── source
├── data
├── outputs
├── shell scripts
├── GWAS results
└── figures
```

- The `data` directory contains all the raw data that was used for this project.
- The `outputs` directory contains the files generated from the raw data through data curation.
- The `source` directory contains the scripts used for data curation and performing various analyses of this project.
- The `shell scripts` directory contains the shell scripts to run code over the command line.
- The `figures` directory contains the intermediary figures generated for this project. Final figures were assembled using Adobe Illustrator. 


### 3c. Description of files

####Source####

1. `phenotype_curation.Rmd` Code for curating phenotype table and outputting a list of phenotypes.
2. `pop_analyses` Code for visualizing prinicipal components analysis.
3. `gwas_script.Rmd` script to create shell scripts to filter phenotype and genotype data, create kinship matrices, run PCA.
4. `simple_m.R script` to run Simple M to calculate the effective number of markers.
5. `mlmm_gwas_batch_script_final.R` code to run the MLMM GWAS.
6. `correlations.R` code to run and visualize correlations between phenotypes of interest.
7. `top_snp_genotypes.Rmd` code to extract the genotypes of the top SNP hits from the GWAS of interest.  
8. `ripening_model.Rmd` Script to create the ripening model figure.
9. `boxplots_manhattans.Rmd` Code for plotting manhattan plots and boxplots.
10. `zoom_plots` Script to plot zoom ins of manhattan plots with gene annotations.
11. `snp_variance` Script to calculate the proportion of variance explained by the top GWAS snps.


####Data####

1. `pheno_meta_data.csv` Phenotype data and meta data for all accessions in the ABC from [Watts et al. 2021](https://nph.onlinelibrary.wiley.com/doi/full/10.1002/ppp3.10211).
2. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001.nosex` File with the apple IDS that have genetic data.
3. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001.frq` Minor allele frequency file for ABC snps.
4. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned1.txt` PC values from PCA, contig snps removed and LD pruned.
5. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_pruned2.txt` Eigen values from PCA, contig snps removed and LD pruned.
6. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_het.het` Plink file with heterozygosity per individual.
7. `top_snps.txt` file with the names of the top SNP hits from the GWAS of interest.
8. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_top_gwas_snps.raw` Genotype file that has been subset to only include the top SNPs of interest from the GWAS.
9. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_top_gwas_snps.ped` PED genotype file that has been subset to only include the top SNPs of interest from the GWAS.
10. `abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_top_gwas_snps.map` MAP genotype file that has been subset to only include the top SNPs of interest from the GWAS.
11. `gene_models_20170612.gff3.gz` Gene model annotations from the GDDH genome version 1.
12. `vineland_snps.txt` The list of SNPs that were genotyped using HRM.

####Outputs####

1. `geno_pheno_meta_data.csv` Filtered phenotype table that includes 1054 apple IDs that have both phenotype data and genotype data.
2. `pheno_list.txt` List with names of phenotypes from [Watts et al. 2021](https://nph.onlinelibrary.wiley.com/doi/full/10.1002/ppp3.10211).
3. `pheno_list_2017` List with names of phenotypes used for GWAS.
3. `contig_snps_remove.txt` List of unanchored SNPs to be removed from SNP table.
4. `pc1_vs_phenos.csv`The R-square and p-values from the correlation of phenotypes with PC1.
5. `pc2_vs_phenos.csv` The R-square and p-values from the correlation of phenotypes with PC2.
6. `phenotype_sample_sizes.txt` Table with samples sizes for each phenotype.
7. `phenotypes4gwas` Folder that contains two files for each phenotype: one with a list of apple IDs for that phenotype that have trait data and the second a file with the trait measurements per apple ID for that phenotype.
8. `simple_m.out` Output from Simple M package that calculates the effective number of markers.
9. `top_snp_genos.csv` File with the genotypes of the top SNP hits from the GWAS.
10. `gene_annotations` Folder containing the files with gene annotations surrounding the top GWAS hits.
11. `snp_variation.csv` R-square values from LMs with top GWAS SNPs and PCs.
12. `top_snps_pheno_pcs.csv` File with genotypes of the top SNPs from GWAS, pheno data and PCs.
13. `ripening_model_summary.csv` Median values for traits measurements across the genotypic classes at NAC18.1


####Shell Scripts####

1. `abc_pca.sh` script to run PCA for the whole SNP set with TASSEL.
2. `genotype_filtering_plink.sh` script that contains PLINK commands to filter the ABC MAP and PED to only containing apple IDs for a particular phenotype, applies a MAF filter of 0.01, and outputs a MAP and PED for each phenotype.
2. `kinship.sh` script that contains the tassel commands to make a kinship matrix for each phenotype.
3. `pca.sh` script commands to run PCA with tassel for each individual phenotype file.
5. `geno_raw.sh` script with commands to recode PED and MAP files into .raw files for the MLMM gwas.
6. `simple_m.sh` script to run simple_m.R
7. `run_mlmm_gwas_batch_script_final.sh` script to run the MLMM GWAS code.


####GWAS results####

1. `mlmm_pvals` Folder containing files for each phenotype with the SNP p-values from the MLMM GWAS.
2. `mlmm_qq` Folder containing qq-plots from each MLMM GWAS.
3. `mlmm_manhattans` Folder containing manhattan plots from each MLMM GWAS.
4. `rss` Folder containing files for each phenotype with the variance explained by the co-factor SNPs at each step of the MLMM GWAS.
5. `standard_pvals` Folder containing files for each phenotype with the SNP p-values from the standard (MLM) GWAS.
5. `standard_qq` Folder containing qq-plots from each standard (MLM) GWAS.
6. `standard_manhattans` Folder containing manhattan plots from each standard (MLM) GWAS.
7. `sbatch_command_mlmm.txt` commands to excute `run_mlmm_gwas_batch_script_final.sh`.









