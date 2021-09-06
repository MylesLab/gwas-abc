setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/pca")

library(tidyverse)

#now going to try bins based on physical distances instead of number of SNPs.
ld_table <-  read_table2("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_ld.ld")

ld_table <- ld_table %>% mutate(snp_dist = BP_B - BP_A)

#plot out to a mega base without dots
jpeg("abc_ld_decay_all_dots.jpeg", width = 500, height =500)
ggplot(ld_table, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, color = "black", size=0.5)+
  geom_hline(yintercept=0.005396065, linetype="dashed")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()

jpeg("abc_ld_decay_all_dots_mega.jpeg", width = 500, height =500)
ggplot(ld_table, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, color = "black", size=0.5)+
  coord_cartesian(xlim = c(0,1000009))+
  geom_hline(yintercept=0.005396065, linetype="dashed")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()

jpeg("abc_ld_decay_all_dots_025.jpeg", width = 500, height =500)
ggplot(ld_table, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, color = "black", size=0.5)+
  coord_cartesian(ylim = c(0, 0.25))+
  geom_hline(yintercept=0.005396065, linetype="dashed")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()

jpeg("abc_ld_decay_all_dots_mega_025.jpeg", width = 500, height =500)
ggplot(ld_table, aes(x=snp_dist, y= R2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, color = "black", size=0.5)+
  coord_cartesian(xlim = c(0,1000009), ylim = c(0, 0.25))+
  geom_hline(yintercept=0.005396065, linetype="dashed")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()
