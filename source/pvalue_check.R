#load packages
library(tidyverse)

#The TASSEL manhattan plots for certain phenotypes  appear to have weird artifacts and unusually high p-values. We look into why this might be in this script by analyzing the MAF, HWE and p-values

#Plot histograms of tassel pvalues
#Correlate MAF with tassel pvalues
#Correlate HWE pvalues with tassel pvalues
#First take the p-values for a SNP and correlate with MAF.

#load list of phenos
pheno_list <- read_csv("pheno_list.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1
r2 = p = c()
for (i in pheno_list) {
  #load tassel gwas p-values
  trait_pvals <- read_delim(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/tassel_gwas/",i,"_",i,"_geno_filtered_+_",i,"_pheno_reformated_+_",i,"_pca_1_stats.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  trait_pvals <- trait_pvals[-1, ]
  trait_pvals <- trait_pvals[, c(2,7)]
  
  #transform to minus log 10 p-values
  trait_pvals <- trait_pvals %>% mutate(log_pval = -log10(p))
  trait_pvals <- trait_pvals %>% select(Marker, log_pval)
  
  
  #plot
  jpeg(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/plots/hist_",i,".jpeg"), width = 11, height = 5)
  hist(trait_pvals$log_pval, breaks=100000, main = "acidity")
  dev.off()
  
  #load maf
  maf_table <- read_delim(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/",i,"_",i,"_geno_filtered.frq"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  maf_table <- maf_table[, c(2,5)]
  
  #join maf and pvals together
  table <- left_join(maf_table, trait_pvals, by = c("SNP"="Marker"))
  
  #plot
  pheno_maf_cor <- ggplot(table, aes(x=MAF, y=log_pval))+
    geom_point(size=3, stroke=0, alpha=0.8)+
    theme_bw()+
    labs(x="MAF", y="p-value")+
    theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
    stat_smooth(method="lm", col="black", se = FALSE)+
    ggtitle(paste("r2 = ", signif(cor.test(table$MAF, table$log_pval, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(table$MAF, table$log_pval, method = "pearson")$p.value, digits = 2)))
  ggsave(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/plots/",i,".jpeg"), plot = pheno_maf_cor)
  
  #extract the p-values and r2 of the correlation between MAF and p-values
  correlation <- cor.test(table$MAF, table$log_pval, method = "pearson")
  r2[i] = signif(correlation$estimate^2)
  p[i] = signif(correlation$p.val)
}

cor_table = cbind(pheno_list, r2, p)
cor_table <- as.data.frame(cor_table)
cor_table$p <- as.numeric(as.character(cor_table$p))
cor_table$r2 <- as.numeric(as.character(cor_table$r2))
cor_table <- cor_table %>% arrange(desc(r2))

write.csv(cor_table,file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/plots/cor_table.csv", quote = F, row.names = F))
          
          
#plot distributions of 
