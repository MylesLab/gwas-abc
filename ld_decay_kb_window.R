#make an LD decay curve that only go out to a certain window. Plot the loess curve and the dots.

setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/pca")

library(tidyverse)

#try bins based on number of SNPs.
ld_table = read.csv("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_ld.ld", sep="", header = T)

ld_table <- ld_table %>% mutate(snp_dist = BP_B - BP_A)

####################################################################################

ld_table_50kb <- ld_table %>% filter(snp_dist < 50000)
#plot out to 50kb
pdf("ld_decay_50kb.pdf", width = 5, height =5)
ggplot(ld_table_50kb, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, size=0.5)+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()


ld_table_10kb <- ld_table %>% filter(snp_dist < 10000)
#plot out to 10kb
pdf("ld_decay_10kb.pdf", width = 5, height =5)
ggplot(ld_table_10kb, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, size=0.5)+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()

ld_table_1kb <- ld_table %>% filter(snp_dist < 1000)
#plot out to 1kb
pdf("ld_decay_1kb.pdf", width = 5, height =5)
ggplot(ld_table_1kb, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, size=0.5)+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()


ld_table_500bp <- ld_table %>% filter(snp_dist < 500)
#plot out to 1kb
pdf("ld_decay_500bp.pdf", width = 5, height =5)
ggplot(ld_table_500bp, aes(x=snp_dist, y= R2))+ 
  geom_point()+
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, size=0.5)+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()

