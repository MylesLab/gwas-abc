#We determined that these weird artifacts can't be explained by MAF or HWE filtering. We though maybe this could be due to the compression of the kinship matrix that takes place in TASSEL. Therefore we tried running tassel again just for phenolic content with no compression (see code in alt_tassel_manhattan.R) to see if it made a difference. Now I will correlate the p-values from the orginal tassel to the no compression tassel p-values.

tassel <- read_table2("original_tpc_tpc_geno_filtered_+_tpc_pheno_reformated_+_tpc_pca_1_stats.txt")
nocomp_tassel <- read_table2("tpc_tpc_geno_filtered_+_tpc_pheno_reformated_+_tpc_pca_1_stats.txt")

tassel <- tassel[-1, ]
nocomp_tassel <- nocomp_tassel[-1, ]

tassel <- tassel %>% select(Marker, p)
nocomp_tassel <- nocomp_tassel %>% select(Marker, p)

tassel$p <- as.numeric(as.character(tassel$p))
nocomp_tassel$p <- as.numeric(as.character(nocomp_tassel$p))

tassel <- tassel %>% mutate(log_pval = -log10(p))
nocomp_tassel <- nocomp_tassel %>% mutate(log_pval = -log10(p))


nocomp_tassel <- nocomp_tassel %>% rename(nocomp_log_pval = log_pval)


pvalues <- left_join(tassel, nocomp_tassel, by = "Marker")

pvalues <- pvalues %>% select(Marker, log_pval, nocomp_log_pval)

pval_cor <- ggplot(pvalues, aes(x=log_pval, y=nocomp_log_pval))+
  geom_point(size=1, stroke=0, alpha=0.8)+
  theme_bw()+
  coord_fixed(xlim=c(0,25), ylim=c(0,25))+
  labs(x="Tassel with compression -log10(p)", y="TASSEL without compression -log10(p)")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  ggtitle(paste("r2 = ", signif(cor.test(pvalues$log_pval, pvalues$nocomp_log_pval, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(pvalues$log_pval, pvalues$nocomp_log_pval, method = "pearson")$p.value, digits = 2)))
ggsave("tassel_pval_cor.jpeg", plot = pval_cor)


