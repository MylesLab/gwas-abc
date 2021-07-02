#set working directory
setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/tassel_gwas/tassel_manhattans")

#load packages
library(tidyverse)
library(qqman)
library(scales)

#Manhattan plot for GWAS in TASSEL
#load list of phenos
pheno_list <- read_csv("outputs/pheno_list.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1

#load full pheno names
full_pheno_names <- read_delim("data/full_pheno_names.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
short <- full_pheno_names$name
long <- full_pheno_names$full_name
names(long) <- short

#load SNPs to highlight
vineland_snps <- read_csv("data/vineland_snps.txt")
vineland_snps <- vineland_snps$snp

#set chromosome names
chromos <- c(1:17, "R")
chromos <- as.character(chromos)

for (i in pheno_list) {
  gwas_output <- read_delim(file = paste0(i,"_",i,"_geno_filtered_+_",i,"_pheno_reformated_+_",i,"_pca_1_stats.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  gwas_output <- as.data.frame(gwas_output)
  gwas_output <- gwas_output[-1, ]
  gwas_output$Chr <- as.numeric(as.character(gwas_output$Chr))
  gwas_output$Pos <- as.numeric(as.character(gwas_output$Pos))
  gwas_output$p <- as.numeric(as.character(gwas_output$p))
  #make chr 0 be chr R
  gwas_output[which(gwas_output$Chr == 0), "Chr"] <- 18
  jpeg(file=paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/tassel_gwas/tassel_manhattans/",i,"_tassel.jpeg"), width = 600, height = 300)
  par(mfcol=c(1,2))
  manhattan_table = gwas_output[,c("Marker","Chr", "Pos", "p")]
  colnames(manhattan_table) = c("SNP","CHR", "BP", "P")
  manhattan_table=manhattan_table[c("BP","CHR", "P","SNP")]
  max_p <- -log10(min(manhattan_table$P)*0.01)
  jpeg(file="firmness_harv_17.jpeg", width = 600, height = 300)
  manhattan(manhattan_table, main = long[i], suggestiveline=F, genomewideline=-log10(0.05/211156), highlight = vineland_snp), chrlabs = chromos, ylim = c(0,max_p))
  dev.off()
  qq(manhattan_table$P, main = long[i])
  dev.off()
}
