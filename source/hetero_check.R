#load packages
library(tidyverse)

#load list of phenos
pheno_list <- read_csv("pheno_list.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1

for (i in pheno_list) {
  
  #load TASSEL gwas p-values
  trait_pvals <- read_delim(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/tassel_gwas/",i,"_",i,"_geno_filtered_+_",i,"_pheno_reformated_+_",i,"_pca_1_stats.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  trait_pvals <- trait_pvals[-1, ]
  trait_pvals <- trait_pvals[, c(2,7)]
  trait_pvals$p <- as.numeric(as.character(trait_pvals$p))
  
  #transform to minus log 10 p-values
  trait_pvals <- trait_pvals %>% mutate(log_pval = -log10(p))
  trait_pvals <- trait_pvals %>% select(Marker, log_pval)
  
  #load maf
  maf_table <- read_table2(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/",i,"_geno_filtered.frq"))
  maf_table <- maf_table[, c(2,5)]
  
  
  #join maf and pvals together
  table <- left_join(maf_table, trait_pvals, by = c("SNP"="Marker"))
  
  
  #load HWE data.
  hardy_p <- read_table2(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/hardy/",i,"_geno_filtered.hwe"))
  hardy_p <- hardy_p %>% mutate(hwe_log_p = -log10(P)) %>% rename(heteros = 'O(HET)')
  
  #make list of SNPs that are more than 90% heterozygous across individuals.
  hetero_list <- hardy_p %>% filter(heteros > 0.9000)
  hetero_list <- hetero_list %>% select(SNP) %>% mutate(snp2 = SNP)
  write.table(hetero_list, file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/heteros/",i,"_heteros.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
  
  hardy_p <- hardy_p %>% select(SNP, heteros, hwe_log_p)
  
  #join HWE table with the big table.
  table <- left_join(table, hardy_p, by = "SNP")
  
  #plot
  het_maf_cor <- ggplot(table, aes(x=heteros, y=MAF))+
    geom_point(size=1, stroke=0, alpha=0.8)+
    theme_bw()+
    labs(x="Proportion of heterozygotes", y="MAF")+
    theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
    ggtitle(paste("r2 = ", signif(cor.test(table$heteros, table$MAF, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(table$heteros, table$MAF, method = "pearson")$p.value, digits = 2)))
  ggsave(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/heteros/cor_plots/het_maf_cor_",i,".jpeg"), plot = het_maf_cor)
  
  het_pval_cor <- ggplot(table, aes(x=heteros, y=log_pval))+
    geom_point(size=1, stroke=0, alpha=0.8)+
    theme_bw()+
    labs(x="Proportion of heterozygotes", y="TASSEL -log10(p)")+
    theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
    ggtitle(paste("r2 = ", signif(cor.test(table$heteros, table$log_pval, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(table$heteros, table$log_pval, method = "pearson")$p.value, digits = 2)))
  ggsave(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/heteros/cor_plots/het_pval_cor_",i,".jpeg"), plot = het_pval_cor)
  
  hwe_maf_cor <- ggplot(table, aes(x=hwe_log_p, y=MAF))+
    geom_point(size=1, stroke=0, alpha=0.8)+
    theme_bw()+
    labs(x="HWE -log10(p)", y="MAF")+
    theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
    ggtitle(paste("r2 = ", signif(cor.test(table$hwe_log_p, table$MAF, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(table$hwe_log_p, table$MAF, method = "pearson")$p.value, digits = 2)))
  ggsave(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/heteros/cor_plots/hwe_maf_cor_",i,".jpeg"), plot = hwe_maf_cor)
  
  

}