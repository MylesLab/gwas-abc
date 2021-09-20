#set working directory
setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas")

#load packages
library(tidyverse)
library(qqman)
library(scales)

#Manhattan plot for GWAS in TASSEL

#load SNPs to highlight
vineland_snps <- read_csv("vineland_snps.txt")
vineland_snps <- vineland_snps$snp

#set chromosome names
chromos <- c(1:17, "R")
chromos <- as.character(chromos)


gwas_output <- read_delim(file = "/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/alt_tassel_gwas/tpc_tpc_geno_filtered_+_tpc_pheno_reformated_+_tpc_pca_1_stats.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
gwas_output <- as.data.frame(gwas_output)
gwas_output <- gwas_output[-1, ]
gwas_output$Chr <- as.numeric(as.character(gwas_output$Chr))
gwas_output$Pos <- as.numeric(as.character(gwas_output$Pos))
gwas_output$p <- as.numeric(as.character(gwas_output$p))
#make chr 0 be chr R
gwas_output[which(gwas_output$Chr == 0), "Chr"] <- 18
manhattan_table = gwas_output[,c("Marker","Chr", "Pos", "p")]
colnames(manhattan_table) = c("SNP","CHR", "BP", "P")
manhattan_table=manhattan_table[c("BP","CHR", "P","SNP")]
max_p <- -log10(min(manhattan_table$P)*0.01)

jpeg(file="/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/alt_tassel_gwas/tpc_alt_manhattan.jpeg", width = 800, height = 400)
manhattan(manhattan_table, main = "Phenolics", suggestiveline=F, genomewideline=-log10(0.05/211156), highlight = vineland_snps, chrlabs = chromos, ylim = c(0,max_p))
dev.off()

jpeg(file="/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/alt_tassel_gwas/tpc_alt_qq.jpeg", width = 500, height = 500)
qq(manhattan_table$P, main = "Phenolics")
dev.off()

