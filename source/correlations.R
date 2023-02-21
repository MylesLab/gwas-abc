library(tidyverse)
library(viridis)
library(GGally)


#make correlation plot of phenoytypes with NAC hit.

geno_pheno_meta_data <- read_csv("../outputs/geno_pheno_meta_data.csv")

abc_phenos <- geno_pheno_meta_data %>% select(date_jul_17_harv, firmness_avg_17_harv, percent_firmness_avg_17, firmness_avg_17_stor, brix_17_harv, percent_brix_17, brix_17_stor, juiciness_16_harv)

#make plot
nac_phenos_correlations <- ggpairs(abc_phenos,
        upper = list(continuous = wrap("cor", size = 2.5)),
        lower = list(continuous = wrap("points",size = 0.1, alpha = 0.5)))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank())
ggsave("../figures/nac_phenos_correlations.pdf", plot = nac_phenos_correlations)

