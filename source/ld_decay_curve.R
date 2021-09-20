setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/pca")

library(tidyverse)

back_ld <- read_table2("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_back_ld.ld")
back_ld <- subset(back_ld, (CHR_A != CHR_B))
avg_back_ld <- mean(back_ld$R2, na.rm = T)

#now going to try bins based on physical distances instead of number of SNPs.
ld_table = read.csv("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_ld.ld", sep="", header = T)
snp_count = length(unique(c(as.character(ld_table[, "SNP_A"]), as.character(ld_table[, "SNP_B"]))))
snp_dist = as.numeric(ld_table[, "BP_B"] - ld_table[, "BP_A"])
r2 = as.numeric(ld_table[, "R2"])
snp_r2 <- cbind(snp_dist, r2)
samp_bins = seq(0, 1000000, by = 1000)
median_r2 = mean_r2 = c()
for (i in 1:(length(samp_bins) - 1)) {
  index = which(snp_dist>samp_bins[i] & snp_dist<=samp_bins[i+1]) #elements where snp dist is within 
  median_r2[i] = median(r2[index])
  mean_r2[i] = mean(r2[index])
}
ld_stats <- cbind(median_r2, mean_r2, samp_bins[2:length(samp_bins)])
ld_stats <- as.data.frame(ld_stats)

#plot out to a mega base without dots.
pdf("abc_ld_decay_025.pdf", width = 5, height =5)
ggplot(ld_stats, aes(x=V3, y= mean_r2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, color = "black", size=0.5)+
  coord_cartesian(xlim = c(0,1000009), ylim = c(0, 0.25))+
  geom_hline(yintercept=avg_back_ld, linetype="dashed")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()


