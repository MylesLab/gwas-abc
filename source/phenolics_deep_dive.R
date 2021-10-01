#look at kendra's phenolic data to see if the overlapping apple accessions from her study and the ABC GWAS study that were easured for phenolic content have similar levels of phenolics. We used a different assay to measure phenolics so we would like to see if they are correlated and therefore if we would expect to see the same GWAS hit for phenolics.

library(readxl)
library(tidyverse)

phenolics_2014 <- read_excel("data/41438_2019_190_MOESM2_ESM (2).xls")
phenolics_2016 <- read_excel("data/41438_2019_190_MOESM3_ESM (4).xls")

#keep total phenolics
phenolics_2014 <- phenolics_2014 %>% select(ID, total_phenolics)
phenolics_2016 <- phenolics_2016 %>% select(ID, total_phenolics)

#load pheno data
geno_pheno_meta_data <- read_csv("outputs/geno_pheno_meta_data.csv")

phenolics_2014 <- left_join(phenolics_2014, geno_pheno_meta_data, by = c("ID" = "PLANTID"))
phenolics_2016 <- left_join(phenolics_2016, geno_pheno_meta_data, by = c("ID" = "PLANTID"))


tpc_2014 <- ggplot(phenolics_2014, aes(x=tpc, y=total_phenolics))+
  geom_point(size=3, stroke=0, alpha=0.8)+
  theme_bw()+
  labs(x="Folin–Ciocalteu assay TPC", y="HPLC TPC 2014")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(phenolics_2014$tpc, phenolics_2014$total_phenolics, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(phenolics_2014$tpc, phenolics_2014$total_phenolics, method = "pearson")$p.value, digits = 2)))
ggsave("figures/tpc_2014.jpeg", plot = tpc_2014)

tpc_2016 <- ggplot(phenolics_2016, aes(x=tpc, y=total_phenolics))+
  geom_point(size=3, stroke=0, alpha=0.8)+
  theme_bw()+
  labs(x="Folin–Ciocalteu assay TPC", y="HPLC TPC 2016")+
  theme(panel.border = element_blank(), axis.text=element_text(colour = "black", size=10), axis.title=element_text(size = 10, face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(size = 9))+
  stat_smooth(method="lm", col="black", se = FALSE)+
  ggtitle(paste("r2 = ", signif(cor.test(phenolics_2016$tpc, phenolics_2016$total_phenolics, method = "pearson")$estimate^2, digits = 3), " p-value = ", signif(cor.test(phenolics_2016$tpc, phenolics_2016$total_phenolics, method = "pearson")$p.value, digits = 2)))
ggsave("figures/tpc_2016.jpeg", plot = tpc_2016)


