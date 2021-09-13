library(tidyverse)

#Some phenotypes in the TASSEL GWAS had numerous p-values clustered together which appeared strange and looked like weird artifacts on the QQ and Manhattan plots. We compared these TASSEL p-values to MAF and HWE p-values in the "pvalue_check.R" script. Here we will continue a couple other checks.

#load in sample sizes for each phenotype we ran GWAS for.
phenotype_sample_sizes <- read_delim("outputs/phenotype_sample_sizes.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

#Load in the p-values from TASSEL

#Bin p-values and then calculate the median