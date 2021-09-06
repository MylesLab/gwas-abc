library(tidyverse)
library(viridis)
library(GGally)


#make correlation of 
geno_pheno_meta_data <- read_csv("outputs/geno_pheno_meta_data.csv")

abc_phenos <- geno_pheno_meta_data %>% select(date_jul_17_harv, firmness_avg_17_harv, percent_firmness_avg_17, firmness_avg_17_stor, percent_brix_17, brix_17_stor)




#make correlation plot of phenoytypes with NAC hit.

geno_pheno_meta_data <- read_csv("outputs/geno_pheno_meta_data.csv")

abc_phenos <- geno_pheno_meta_data %>% select(time_ripen_2017, date_jul_17_harv, firmness_avg_17_harv, percent_brix_17, brix_17_stor, juiciness_16_harv)


abc_phenos[1:6] <- lapply(abc_phenos[1:6], as.character)
abc_phenos[1:6] <- lapply(abc_phenos[1:6], as.numeric)

#change from tbl to base dataframe.
abc_phenos <- as.data.frame(abc_phenos)

#create an empty matrix for correlations, make matrix with no data and just col and rows from final pheno table.
pairwise_pheno_correlations=matrix(,ncol(abc_phenos), ncol(abc_phenos))
rownames(pairwise_pheno_correlations)=colnames(abc_phenos)
colnames(pairwise_pheno_correlations)=colnames(abc_phenos)
#marix for pvalues.
pairwise_pheno_correlations_pval=matrix(,ncol(abc_phenos), ncol(abc_phenos))
rownames(pairwise_pheno_correlations_pval)=colnames(abc_phenos)
colnames(pairwise_pheno_correlations_pval)=colnames(abc_phenos)

#Treat the data as all quantitative data and run pearson's correlation.
for (i in 1:ncol(abc_phenos)) {
  phenoname_x = colnames(abc_phenos)[i]
  for (j in 1:ncol(abc_phenos)) {
    phenoname_y = colnames(abc_phenos)[j]
    pairwise_pheno_correlations[j,i]=cor.test(abc_phenos[,i], abc_phenos[,j], method = "pearson")$estimate^2
    pairwise_pheno_correlations_pval[j,i]= cor.test(abc_phenos[,i], abc_phenos[,j], method = "pearson")$p.value
  }
}

#bonneferonni correct
pairwise_pheno_correlations_pval[upper.tri(pairwise_pheno_correlations_pval)] = p.adjust(pairwise_pheno_correlations_pval[upper.tri(pairwise_pheno_correlations_pval)], method = "bonferroni")
pairwise_pheno_correlations_pval[lower.tri(pairwise_pheno_correlations_pval)] = p.adjust(pairwise_pheno_correlations_pval[lower.tri(pairwise_pheno_correlations_pval)], method = "bonferroni")

nac_phenos_correlations <- ggpairs(abc_phenos,
        upper = list(continuous = wrap("cor", size = 2.5)),
        lower = list(continuous = wrap("points",size = 0.1, alpha = 0.5)))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank())
ggsave("figures/nac_phenos_correlations.pdf", plot = nac_phenos_correlations)


##Soluble solids

ssc_firmness <- ggplot(geno_pheno_meta_data, aes(x=brix_17_stor, y=firmness_avg_17_harv))+
  geom_point(aes(colour=percent_brix_17),size=3, stroke=0, alpha=0.8)+
  theme_bw()+
  coord_fixed()+
  scale_color_viridis(name = "Change in SSC (%)", option = "magma", direction = -1)+
  labs(x="SSC after storage", y="Firmness at harvest")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(geno_pheno_meta_data$firmness_avg_17_harv, geno_pheno_meta_data$brix_17_stor, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(geno_pheno_meta_data$firmness_avg_17_harv, geno_pheno_meta_data$brix_17_stor, method = "pearson")$p.value, digits = 2)))
ggsave("figures/ssc_firmness.jpeg", plot = ssc_firmness)

change_ssc_firmness <- ggplot(geno_pheno_meta_data, aes(x=percent_brix_17, y=firmness_avg_17_harv))+
  geom_point(aes(colour=brix_17_stor),size=3, stroke=0, alpha=0.8)+
  theme_bw()+
  scale_color_viridis(name = "SSC after storage", option = "magma", direction = -1)+
  labs(x="Change in SSC during storage", y="Firmness at harvest")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(geno_pheno_meta_data$firmness_avg_17_harv, geno_pheno_meta_data$percent_brix_17, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(geno_pheno_meta_data$firmness_avg_17_harv, geno_pheno_meta_data$percent_brix_17, method = "pearson")$p.value, digits = 2)))
ggsave("figures/change_ssc_firmness.jpeg", plot = change_ssc_firmness)

ssc_harvest_firmness <- ggplot(geno_pheno_meta_data, aes(x=brix_17_harv, y=firmness_avg_17_harv))+
  geom_point(size=3, stroke=0, alpha=0.6)+
  theme_bw()+
  coord_fixed()+
  labs(x="SSC at harvest", y="Firmness at harvest")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(geno_pheno_meta_data$firmness_avg_17_harv, geno_pheno_meta_data$brix_17_harv, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(geno_pheno_meta_data$firmness_avg_17_harv, geno_pheno_meta_data$brix_17_harv, method = "pearson")$p.value, digits = 2)))
ggsave("figures/ssc_harvest_firmness.jpeg", plot = ssc_harvest_firmness)

##SSC and harvest date

ssc_date <- ggplot(geno_pheno_meta_data, aes(x=date_jul_17_harv, y=brix_17_stor))+
  geom_point(aes(colour=percent_brix_17),size=3, stroke=0, alpha=0.8)+
  theme_bw()+
  scale_color_viridis(name = "Change in SSC (%)", option = "magma", direction = -1)+
  labs(x="Harvest date", y="SSC after storage")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(geno_pheno_meta_data$brix_17_stor, geno_pheno_meta_data$date_jul_17_harv, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(geno_pheno_meta_data$brix_17_stor, geno_pheno_meta_data$date_jul_17_harv, method = "pearson")$p.value, digits = 2)))
ggsave("figures/ssc_date.jpeg", plot = ssc_date)

change_ssc_date <- ggplot(geno_pheno_meta_data, aes(x=date_jul_17_harv, y=percent_brix_17))+
  geom_point(aes(colour=brix_17_stor),size=3, stroke=0, alpha=0.8)+
  theme_bw()+
  scale_color_viridis(name = "SSC after storage", option = "magma", direction = -1)+
  labs(x="Harvest date", y="Change in SSC during storage")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(geno_pheno_meta_data$percent_brix_17, geno_pheno_meta_data$date_jul_17_harv, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(geno_pheno_meta_data$percent_brix_17, geno_pheno_meta_data$date_jul_17_harv, method = "pearson")$p.value, digits = 2)))
ggsave("figures/change_ssc_date.jpeg", plot = change_ssc_date)


###Juiciness

juicy_harvest <- ggplot(geno_pheno_meta_data, aes(x=date_jul_17_harv, y=juiciness_16_harv))+
  geom_point(,size=3, stroke=0, alpha=0.6)+
  theme_bw()+
  labs(x="Harvest date (Julian day)", y="Juiciness")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(geno_pheno_meta_data$date_jul_17_harv, geno_pheno_meta_data$juiciness_16_harv, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(geno_pheno_meta_data$date_jul_17_harv, geno_pheno_meta_data$juiciness_16_harv, method = "pearson")$p.value, digits = 2)))
ggsave("figures/juicy_harvest.jpeg", plot = juicy_harvest)
