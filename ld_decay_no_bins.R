#LD decay with no bins

setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/pca")

library(tidyverse)

#try bins based on number of SNPs.
ld_table = read.csv("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_ld.ld", sep="", header = T)

ld_table <- ld_table %>% mutate(snp_dist = BP_B - BP_A)

#plot
pdf("ld_decay_no_bins.pdf", width = 5, height =5)
ggplot(ld_table, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, size=0.5)+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()
