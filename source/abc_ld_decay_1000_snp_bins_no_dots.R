setwd("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/pca")

library(tidyverse)

#try bins based on number of SNPs.
ld_table = read.csv("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_ld.ld", sep="", header = T)

#loop through 100 rows at a time and calculate the median and mean R square
ld_table = as.matrix(ld_table)
snp_dist = as.numeric(as.numeric(ld_table[, "BP_B"]) - as.numeric(ld_table[, "BP_A"]))
r2 = as.numeric(ld_table[, "R2"])
r2 = r2[order(snp_dist)]
snp_dist = snp_dist[order(snp_dist)]
samp_bins = seq(1, length(snp_dist), by = 1000)
median_r2 = c()
mean_dist = c()
for (i in 1:(length(samp_bins) - 1)) {
  mean_dist[i] = mean(snp_dist[samp_bins[i]:(samp_bins[i + 1] - 1)])
  median_r2[i] = median(r2[samp_bins[i]:(samp_bins[i + 1] - 1)])
}

table = cbind(mean_dist, median_r2)

table <- as.data.frame(table)

#plot
pdf("abc_ld_decay_1000_snp_bins_no_dots.pdf", width = 5, height =5)
ggplot(table, aes(x=mean_dist, y= median_r2))+ 
  geom_smooth(colour = "darkgrey", method = "loess", se=FALSE, size=0.5)+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  labs(x="Inter-SNP distance (bp)", y="LD - r squared")+
  theme_classic()
dev.off()



