
library(tidyverse)

back_ld <- read_table2("abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001_noctgs_back_ld.ld")
back_ld <- subset(back_ld, (CHR_A != CHR_B))
avg_back_ld <- mean(back_ld$R2, na.rm = T)

print(avg_back_ld)