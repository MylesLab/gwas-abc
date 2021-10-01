#I am troubleshooting the MLMM script. I will try making a pheno list file with just one phenotype to see if if the code can run through at least one phenotype. I made two files tpc.txt and  with just "tpc" listed as a phenotype.

#set working directory
setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas")

#load packages
library(tidyverse)
library(emma)
library(mlmm)
library(qqman)
library(scales)

###################################################################################################
#load list of phenos
pheno_list <- read_csv("tpc.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1

#load full pheno names
full_pheno_names <- read_delim("tpc_full_name.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
short <- full_pheno_names$name
long <- full_pheno_names$full_name
names(long) <- short

print(paste0("full_pheno_names cols=",ncol(full_pheno_names)))
print(paste0("full_pheno_names rows=",nrow(full_pheno_names)))

###################################################################################################

for (i in pheno_list) {
  
  #load phenotype data
  pheno_data <- read_delim(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/phenotype_data/", i,"_phenotype.txt"), 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)
  pheno <- pheno_data$X2
  names(pheno) <- pheno_data$X1
  
  print(paste0("pheno_data cols=",ncol(pheno_data)))
  print(paste0("pheno_data rows=",nrow(pheno_data)))
  
  #load genetic data
  geno_dat <- read.table(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/geno_raw/",i,"_geno_filtered.raw"), header=T)
  rownames(geno_dat) <- geno_dat[,1]
  geno_dat <- geno_dat[,7:ncol(geno_dat)]
  
  print(paste0("geno_dat cols=",ncol(geno_dat)))
  print(paste0("geno_dat rows=",nrow(geno_dat)))
  
  #save SNP names from raw file to be able to rename the map file with the same names.
  snp_names <- colnames(geno_dat)
  #get rid of the 'X' at the beginning of each SNP name.
  snp_names <- as.data.frame(snp_names)
  snp_names <- snp_names %>%
    rename(new_snp = snp_names) %>%
    mutate_all(~as.character(.)) %>%
    mutate(trimmed = str_sub(new_snp, start = 2, end= -3))
  
  print(paste0("snp_names cols=",ncol(snp_names)))
  print(paste0("snp_names rows=",nrow(snp_names)))
  
  #make it a matrix
  geno_dat <- as.matrix(geno_dat)
  cat(i," number of SNPs: ", ncol(geno_dat), "\n")
  cat(i," number of apple IDs: ", nrow(geno_dat), "\n")
  
  #load kinship
  kinship <-read_delim(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/kinship/",i,"_kinship.txt"), 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE, skip = 3)
  kinship <- as.data.frame(kinship)
  rownames(kinship) <- kinship[,1]
  kinship <- kinship[,2:ncol(kinship)]
  colnames(kinship) <- rownames(kinship)
  kinship <- as.matrix(kinship)
  
  print(paste0("kinship cols=",ncol(kinship)))
  print(paste0("kinship rows=",nrow(kinship)))
  
  #load PCs
  pcs <- read_delim(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/pca/",i,"_pca_1.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE, 
                    skip = 2)
  pcs <- as.data.frame(pcs)
  rownames(pcs) <- pcs[,1]
  pcs <- pcs[,2:ncol(pcs)]
  colnames(pcs)[c(1:5)] <- c(1:5)
  pcs <- as.matrix(pcs)
  
  print(paste0("pcs cols=",ncol(pcs)))
  print(paste0("pcs rows=",nrow(pcs)))
  
  #map file
  map <- read_delim(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/",i,"_geno_filtered.map"), 
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
  
  map <- map %>%
    select(X1, X2, X4) %>%
    rename(SNP=X2, Chr=X1, Pos=X4) %>%
    select(SNP, Chr, Pos)
  
  #merge with snp names from raw file
  map <- left_join(snp_names, map, by = c("trimmed"="SNP"))
  
  map <- map %>%
    rename(SNP = new_snp) %>%
    mutate(Chr = replace(Chr, Chr == 0, 18)) %>%
    select(SNP, Chr, Pos)
  
  map <- as.data.frame(map)
  
  print(paste0("map rows=",nrow(map)))
  print(paste0("map cols=",ncol(map)))
  
  #run gwas
  gwas_output <- mlmm_cof(Y=pheno, X=geno_dat, cofs=pcs, K=kinship,
                          nbchunks=2, maxsteps=10)
  
  #plot Manhattan for standard mlm
  jpeg(file=paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/mlmm_gwas/standard_manhattans/",i,"_standard_mlm.jpeg"), width = 800, height = 400)
  plot_fwd_GWAS(gwas_output,1,map,1,main=paste0(long[i]," (standard MLM)"), abline(h=(-log10(0.05/211156)), col="black", lty=2))
  dev.off()
  
  #plot Manhattan for optimal mlmm step
  jpeg(file=paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/mlmm_gwas/mlmm_manhattans/",i,"_mlmm.jpeg"), width = 800, height = 400)
  plot_opt_GWAS(gwas_output,'extBIC',map,1,main=paste0(long[i]," (optimal MLMM)"), abline(h=(-log10(0.05/211156)), col="black", lty=2))
  dev.off()
  
  #p-values from standard mlm
  pvals_std <- gwas_output[["pval_step"]][[1]][["out"]]
  colnames(pvals_std)[2] <- "p-value"
  write.table(pvals_std, file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/mlmm_gwas/standard_pvals/",i,"_standard_pvals.csv"), quote = F, sep = ",", row.names = F)
  
  #p-values from optimal mlmm
  pvals_mlmm <- gwas_output[["opt_extBIC"]][["out"]]
  colnames(pvals_mlmm)[2] <- "p-value"
  write.table(pvals_mlmm, file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/mlmm_gwas/mlmm_pvals/",i,"_mlmm_pvals.csv"), quote = F, sep = ",", row.names = F)
  
  #plot QQ for standard mlm
  jpeg(file=paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/mlmm_gwas/standard_qq/",i,"_standard_qq.jpeg"), width = 500, height = 500)
  qq(pvals_std[,2], main = paste0(long[i]," QQ plot (standard MLM)"))
  dev.off()
  
  #plot QQ for mlmm
  jpeg(file=paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/mlmm_gwas/mlmm_qq/",i,"_mlmm_qq.jpeg"), width = 500, height = 500)
  qq(pvals_mlmm[,2], main = paste0(long[i]," QQ plot (optimal MLMM)"))
  dev.off()
  
  #save RSS output
  rss <- as.data.frame(gwas_output[["RSSout"]])
  write.table(rss, file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/mlmm_gwas/rss/",i,"_rss.csv"), quote = F, sep = ",", row.names = F)
}
