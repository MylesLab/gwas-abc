library(tidyverse)

#Some phenotypes in the TASSEL GWAS had numerous p-values clustered together which appeared strange and looked like weird artifacts on the QQ and Manhattan plots. We compared these TASSEL p-values to MAF and HWE p-values in the "pvalue_check.R" script. Here we will continue a couple other checks.

#Load in the p-values from TASSEL

#Bin p-values and then calculate the count, the median MAF and the median HWE for each bin.

#load hardy p-values.

hardy_p <- read_table2("data/acidity_16_harv_geno_filtered.hwe")
