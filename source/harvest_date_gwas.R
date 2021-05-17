setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/harvest_prelim")

library(tidyverse)
library(emma)
library(mlmm)
library(qqman)
library(scales)

#KINSHIP
kinship <- read_delim("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001_kinship.txt", 
              "\t", escape_double = FALSE, col_names = FALSE, 
              trim_ws = TRUE, skip = 3)
kinship <- as.data.frame(kinship)
order <- kinship$X1
rownames(kinship) <- kinship[,1]
kinship <- kinship[,2:ncol(kinship)]
colnames(kinship) <- rownames(kinship)
kinship <- as.matrix(kinship)


#PHENOTYPE
phenotype <- read_delim("harvest_date_pheno.txt", 
           "\t", escape_double = FALSE, col_names = FALSE, 
           trim_ws = TRUE)
phenotype <- phenotype[match(order, phenotype$X1),]
pheno_vec <- phenotype$X2
names(pheno_vec) <- phenotype$X1

#GENOTYPE
geno_dat <- read.table("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001.raw", header=T)
geno_dat <- geno_dat[match(order, geno_dat$FID),]
rownames(geno_dat) <- geno_dat[,1]
geno_dat <- geno_dat[,7:ncol(geno_dat)]
#save colnames from genotype table.
snp_names <- colnames(geno_dat)
geno_dat <- as.matrix(geno_dat)

#arrange snp names df
snp_names_tab <- as.data.frame(snp_names)
snp_names_tab <- snp_names_tab %>%
  rename(new_snp = snp_names) %>%
  mutate_all(~as.character(.)) %>%
  mutate(trimmed = str_sub(new_snp, 2, -3))

#PCA
pcs <- read_delim("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001_1.txt", 
              "\t", escape_double = FALSE, trim_ws = TRUE, 
              skip = 2)
pcs <- as.data.frame(pcs)
pcs <- pcs[match(order, pcs$Taxa),]
rownames(pcs) <- pcs[,1]
pcs <- pcs[,2:ncol(pcs)]
pcs <- as.matrix(pcs)

#load snp info
map <- read_delim("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_harvest_date_maf001.map", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
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

pdf(file = "harv_date_std_manhat.pdf")
plot_fwd_GWAS(gwas_output, 1, map,1,main="Harvest Date (standard MLM)", abline(h=(-log10(0.05/255664)), col="black", lty=2))
dev.off()

std_pvals <- gwas_output[["pval_step"]][[1]][["out"]]
write.table(std_pvals, file = "harv_date_std_pvals.csv", quote = F, sep = ",", row.names = F)

pdf(file = "harv_date_std_qq.pdf")
qqplot <- qq(std_pvals[,2], main = "Harvest Date QQ plot (standard MLM)")
dev.off()



pdf(file = "harv_date_mlmm_manhat.pdf")
plot_opt_GWAS(gwas_output, 'extBIC', map, 1, main="Harvest Date (optimal MLMM)", abline(h=(-log10(0.05/255664)), col="black", lty=2))
dev.off()

mlmm_pvals <- gwas_output[["opt_extBIC"]][["out"]]
write.table(mlmm_pvals, file = "harv_date_mlmm_pvals.csv", quote = F, sep = ",", row.names = F)

pdf(file = "harv_date_mlmm_qq.pdf")
qqplot <- qq(mlmm_pvals[,2], main = "Harvest Date QQ plot (optimal MLMM)")
dev.off()


var <- as.data.frame(gwas_output[["RSSout"]])
write.table(var, file = "harv_date_snp_variation.csv", quote = F, sep = ",", row.names = F)

