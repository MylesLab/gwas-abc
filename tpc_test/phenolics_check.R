library(qqman)
library(scales)


snp_names <- read_csv("tpc_test/tpc_standard_pvals.csv")



snp_names <- as.data.frame(snp_names)
snp_names <- snp_names %>%
  mutate_all(~as.character(.)) %>%
  mutate(trimmed = str_sub(SNP, start = 2, end= -3))

snp_names <- separate(data = snp_names, col = trimmed, into = c("chromo", "bp", "caller"))

snp_names$`p-value` <- as.numeric(as.character(snp_names$`p-value`))

snp_names <- snp_names %>% mutate(log_p = -log10(`p-value`))


mlmm <- read_csv("tpc_test/tpc_mlmm_pvals.csv")

mlmm <- as.data.frame(mlmm)
mlmm <- mlmm %>%
  mutate_all(~as.character(.)) %>%
  mutate(trimmed = str_sub(SNP, start = 2, end= -3))

mlmm <- separate(data = mlmm, col = trimmed, into = c("chromo", "bp", "caller"))

mlmm$`p-value` <- as.numeric(as.character(mlmm$`p-value`))

mlmm <- mlmm %>% mutate(log_p = -log10(`p-value`))


snp_names <- snp_names %>% select(bp, chromo, log_p, caller)


colnames(snp_names) = c("BP","CHR", "P", "SNP")
max_p <- -log10(min(snp_names$P)*0.01)

snp_names$CHR <- as.numeric(as.character(snp_names$CHR))
snp_names$BP <- as.numeric(as.character(snp_names$BP))
snp_names$P <- as.numeric(as.character(snp_names$P))

manhattan(snp_names, main = "Phenolics", suggestiveline=F, genomewideline=-log10(0.05/211156), ylim = c(0,max_p))
