#load packages
library(tidyverse)

#The TASSEL manhattan plots for certain phenotypes appear to have weird artifacts and unusually high p-values. We look into why this might be in this script by analyzing the MAF, HWE and p-values.

#Plot histograms of tassel pvalues
#Correlate MAF with tassel pvalues
#Correlate HWE pvalues with tassel pvalues

#load list of phenos
pheno_list <- read_csv("pheno_list.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1

for (i in pheno_list) {
  
  #load TASSEL gwas p-values
  trait_pvals <- read_delim(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/tassel_gwas/",i,"_",i,"_geno_filtered_+_",i,"_pheno_reformated_+_",i,"_pca_1_stats.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  trait_pvals <- trait_pvals[-1, ]
  trait_pvals <- trait_pvals[, c(2,7)]
  
  #transform to minus log 10 p-values
  trait_pvals <- trait_pvals %>% mutate(log_pval = -log10(p))
  trait_pvals <- trait_pvals %>% select(Marker, log_pval)
  
  #plot histogram of TASSEL p-values
  jpeg(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/plots/pval_hists/tassle_p_hist_",i,".jpeg"), width = 450, height = 500)
  hist(trait_pvals$log_pval, breaks=100000, xlab = "TASSEL -log10(p)", ylab = "Count")
  dev.off()
  
  #load maf
  maf_table <- read_table2(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/",i,"_geno_filtered.frq"))
  maf_table <- maf_table[, c(2,5)]
  
  #plot histogram of MAF values
  jpeg(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/plots/maf_hists/maf_hist_",i,".jpeg"), width = 450, height = 500)
  hist(maf_table$MAF, breaks=100000, xlab = "MAF", ylab = "Count")
  dev.off()
  
  #join maf and pvals together
  table <- left_join(maf_table, trait_pvals, by = c("SNP"="Marker"))

  #plot
  maf_cor <- ggplot(table, aes(x=MAF, y=log_pval))+
    geom_point(size=1, stroke=0, alpha=0.8)+
    theme_bw()+
    labs(x="MAF", y="TASSEL -log10(p)")+
    theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
    ggtitle(paste("r2 = ", signif(cor.test(table$MAF, table$log_pval, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(table$MAF, table$log_pval, method = "pearson")$p.value, digits = 2)))
  ggsave(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/plots/maf_p_cor/maf_cor_",i,".jpeg"), plot = maf_cor)
  
  #load HWE data.
  hardy_dat <- read_table2(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/hardy/",i,"_geno_filtered.hwe"))
  hardy_dat <- hardy_dat %>% mutate(hwe_log_p = -log10(P)) %>% select(SNP, hwe_log_p) 

  jpeg(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/plots/hardy_hists/hardy_hist_",i,".jpeg"), width = 450, height = 500)
  hist(hardy_dat$hwe_log_p, breaks=100, xlab = "HWE -log10(p)", ylab = "Count")
  dev.off()
  
  #join HWE table with the big table.
  table <- left_join(table, hardy_dat, by = "SNP")
    
  #plot
  hardy_cor <- ggplot(table, aes(x=hwe_log_p, y=log_pval))+
  geom_point(size=1, stroke=0, alpha=0.8)+
  theme_bw()+
  labs(x="HWE p-value", y="TASSEL p-value")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  ggtitle(paste("r2 = ", signif(cor.test(table$hwe_log_p, table$log_pval, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(table$hwe_log_p, table$log_pval, method = "pearson")$p.value, digits = 2)))
ggsave(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/plots/hardy_p_cor/hardy_cor_",i,".jpeg"), plot = hardy_cor)
  
}




      