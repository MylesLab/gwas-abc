#Script to reformat phenotype files into the format required by TASSEL.

load(tidyverse)



pheno <- read_csv("outputs/phenotypes4gwas/time_ripen_2017_phenotype.txt", 
         col_names = FALSE)
pheno_header <- read_csv("pheno_header.txt")

colnames(pheno)[1] <- "<Phenotype>"

test <-rbind(pheno_header, pheno, stringAsfactors = T)

write.table(test, "time.txt", sep="\t", row.names = F, quote = F)
