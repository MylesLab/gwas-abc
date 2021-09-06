#set working directory
setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas")

#load packages
library(tidyverse)
library(qqman)
library(scales)

#Manhattan plot for GWAS in TASSEL
#load list of phenos
pheno_list <- read_csv("pheno_list.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1

for (i in pheno_list) {
  gwas_output <- read_delim(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/tassel_gwas/",i,"_",i,"_geno_filtered_+_",i,"_pheno_reformated_+_",i,"_pca_1_stats.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  gwas_output <- as.data.frame(gwas_output)
  gwas_output <- gwas_output[-1, ]
  gwas_output <- gwas_output[, c(2:4,7)]
  gwas_output <- gwas_output[order(gwas_output$p),]
  write.table(gwas_output, file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/tassel_gwas/tassel_pvalues/",i, "_sorted_pvals.txt"), sep = "\t",
              row.names = FALSE)

}
