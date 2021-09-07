#Look at the TASSEL manhattan plots and for phenotypes that appear to have weird artifacts and unusually high p-values, plot the MAF distribution of their genotype files.

#take the p-values for a SNP and correlate with MAF.

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
  #load maf
  maf_table <- read_delim(file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/",i,"_",i,"_geno_filtered.frq"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  maf_table <- maf_table[, c(2,5)]
  #join maf and pvals together
  table <- left_join(maf_table, trait_pvals, by = c("SNP"="Marker"))
  colnames(table)[3] <- "pvalue"
  #plot
  pheno_maf_cor <- ggplot(table, aes(x=MAF, y=pvalue))+
    geom_point(size=3, stroke=0, alpha=0.8)+
    theme_bw()+
    labs(x="MAF", y="p-value")+
    theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
    stat_smooth(method="lm", col="black", se = FALSE)+
    ggtitle(paste("r2 = ", signif(cor.test(table$MAF, table$pvalue, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(table$MAF, table$pvalue, method = "pearson")$p.value, digits = 2)))
  ggsave(paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/plots/",i,".jpeg"), plot = pheno_maf_cor)
  correlation <- cor.test(table$MAF, table$pvalue, method = "pearson")
  r2[i] = signif(correlation$estimate^2)
  p[i] = signif(correlation$p.val)
}

cor_table = cbind(pheno_list, r2, p)
cor_table <- as.data.frame(cor_table)
cor_table$p <- as.numeric(as.character(cor_table$p))
cor_table$r2 <- as.numeric(as.character(cor_table$r2))
cor_table <- cor_table %>% arrange(desc(r2))

write.csv(cor_table,file = paste0("/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/plots/cor_table.csv", quote = F, row.names = F)
