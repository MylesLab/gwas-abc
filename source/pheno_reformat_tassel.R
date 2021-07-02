#Script to reformat phenotype files into the format required by TASSEL.
setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas")

library(tidyverse)

#load phenotype names
pheno_list <- read_csv("pheno_list.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1

#load header
pheno_header <- read_csv("pheno_header.txt")


for (i in pheno_list) {
  pheno <- read_csv(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/phenotype_data/",i,"_phenotype.txt"), 
         col_names = FALSE)

  colnames(pheno)[1] <- "<Phenotype>"

  pheno_reformat <-rbind(pheno_header, pheno)

  pheno_reformat <- pheno_reformat %>% add_row('<Phenotype>' = paste("Taxa",i, sep = "\t"), .before = 2)

  write.table(pheno_reformat, file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/phenotype_data/tassel_format/",i,"_pheno_reformated.txt"), sep="\t", row.names = F, quote = F)
}
