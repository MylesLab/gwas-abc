setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas")

library(tidyverse)
library(emma)
library(mlmm)
library(qqman)
library(scales)

#KINSHIP
kinship <- read_delim("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/kinship/tpc_kinship.txt", 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE, skip = 3)

kinship <- as.data.frame(kinship)
order <- kinship$X1
rownames(kinship) <- kinship[,1]
kinship <- kinship[,2:ncol(kinship)]
colnames(kinship) <- rownames(kinship)
kinship <- as.matrix(kinship)


#PHENOTYPE
phenotype <- read_delim("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/phenotype_data/tpc_phenotype.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

phenotype <- phenotype[match(order, phenotype$X1),]
pheno_vec <- phenotype$X2
names(pheno_vec) <- phenotype$X1

#GENOTYPE
geno_dat <- read.table("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/geno_raw/tpc_geno_filtered.raw", header=T)
  
geno_dat <- geno_dat[match(order, geno_dat$FID),]
rownames(geno_dat) <- geno_dat[,1]
geno_dat <- geno_dat[,7:ncol(geno_dat)]
#save colnames from genotype table.
snp_names <- colnames(geno_dat)
geno_dat <- as.matrix(geno_dat)

#Need to save the snp names from the raw file that have the allele appended to the end and then rename the snps int he map file to match the ones in the raw file. First arrange snp names into df.
snp_names_tab <- as.data.frame(snp_names)
snp_names_tab <- snp_names_tab %>%
  rename(new_snp = snp_names) %>%
  mutate_all(~as.character(.)) %>%
  mutate(trimmed = str_sub(new_snp, 2, -3))

#PCA
pcs <- read_delim("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/pca/tpc_pca_1.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE, 
                   skip = 2)
pcs <- as.data.frame(pcs)
pcs <- pcs[match(order, pcs$Taxa),]
rownames(pcs) <- pcs[,1]
pcs <- pcs[,2:ncol(pcs)]
pcs <- as.matrix(pcs)

#load snp info
map <- read_delim("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/tpc_geno_filtered.map", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

map <- map %>%
  select(X1, X2, X4) %>%
  rename(SNP=X2, Chr=X1, Pos=X4) %>%
  mutate(Chr = replace(Chr, Chr == 0, 18)) %>%
  select(SNP, Chr, Pos)

map <- left_join(snp_names_tab, map, by = c("trimmed"="SNP"))

map <- map %>%
  rename(SNP = new_snp) %>%
  select(SNP, Chr, Pos)

map <- as.data.frame(map)

gwas_output <- mlmm_cof(Y=pheno_vec, X=geno_dat, cofs= pcs, K=kinship,
                        nbchunks=10, maxsteps=10)

pdf(file = "tpc_std_manhat.pdf")
plot_fwd_GWAS(gwas_output, 1, map,1,main="Phenolics (standard MLM)", abline(h=(-log10(0.05/211156)), col="black", lty=2))
dev.off()

std_pvals <- gwas_output[["pval_step"]][[1]][["out"]]
write.table(std_pvals, file = "tpc_std_pvals.csv", quote = F, sep = ",", row.names = F)

pdf(file = "tpc_std_qq.pdf")
qqplot <- qq(std_pvals[,2], main = "Phenolics QQ plot (standard MLM)")
dev.off()



pdf(file = "tpc_mlmm_manhat.pdf")
plot_opt_GWAS(gwas_output, 'extBIC', map, 1, main="Phenolics (optimal MLMM)", abline(h=(-log10(0.05/211156)), col="black", lty=2))
dev.off()

mlmm_pvals <- gwas_output[["opt_extBIC"]][["out"]]
write.table(mlmm_pvals, file = "tpc_mlmm_pvals.csv", quote = F, sep = ",", row.names = F)

pdf(file = "tpc_mlmm_qq.pdf")
qqplot <- qq(mlmm_pvals[,2], main = "Phenolics QQ plot (optimal MLMM)")
dev.off()


var <- as.data.frame(gwas_output[["RSSout"]])
write.table(var, file = "tpc_snp_variation.csv", quote = F, sep = ",", row.names = F)
