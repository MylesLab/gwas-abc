#Figure S1 SNP Stats 

library(tidyverse)
library(patchwork)

snp_freq <- read_table2("../data/abc_combined_maf001_sort_vineland_imputed_pheno_hetero90_maf001.frq")

####SNP DISTRIBUTION####
sort(table(snp_freq$CHR))

#Fewest SNPs are on Chromosome 6 (11026), max number on Chromosome 5 (19613). There are 7757 SNPs on unassembled contigs.

#Look at distributions of SNPs
snp_chromo <- snp_freq %>% mutate(CHR=as.factor(str_replace(CHR, "^0", "R"))) %>% ggplot(aes(x=factor(CHR, level=c("1","2","3", "4","5","6","7","8", "9", "10", "11","12", "13", "14","15","16", "17", "R")))) +geom_bar() + theme_bw() + labs(x="Chromosome", y="Number of SNPs")+theme(axis.title = element_text(face="bold", size=10, colour="black"), axis.text = element_text(colour="black", size=8))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#For each chromosome, what is the average inter-SNP distance 

#Remove unassembled chromosome before calculating 
snp_freq_chr <- snp_freq %>% filter(CHR !="0") %>% separate(SNP, c("chr", "snp", "snp_caller"), sep="_") %>% select(-chr)

inter_snp_dist <- c()
for(i in 1:17){
  snps <- snp_freq_chr %>% filter(CHR==i) %>% arrange(as.numeric(as.character(snp)))  %>% select(snp)
  a <- snps[1:nrow(snps)-1,] %>% pull(snp)
  b <- snps[2:nrow(snps),] %>% pull(snp)
  snp_dist <- as.numeric(as.character(a))-as.numeric(as.character(b)) 
  snp_dist <- abs(snp_dist)
  inter_snp_dist <- append(inter_snp_dist,snp_dist)
}

min(inter_snp_dist)
#1
max(inter_snp_dist)
#1468970
mean(inter_snp_dist)
#2598.298

inter_snp_dist <- inter_snp_dist %>% as.data.frame() %>% rename(snp_dist=".") 

sort(table(inter_snp_dist), dec=T)[1:10]
#1     2     3     4     6     5     7     8     9    10 
#14446  9609  8493  7304  6796  6661  5902  5533  5401  4805 

#There are over 14,000 SNPs only 1 BP apart

#Bin anything greater than 10KB apart 

inter_snp_dist <- inter_snp_dist %>% mutate(snp_dist_bin=ifelse(snp_dist>=10000, 10000,snp_dist))

plot_a <- inter_snp_dist  %>% ggplot(aes(snp_dist_bin)) + geom_histogram(bins=20) + theme_bw()+ labs(x="Inter-SNP Distance (bp)", y="Count")+theme(axis.title = element_text(face="bold", size=10, colour="black"), axis.text = element_text(colour="black", size=8))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_vline(xintercept=2598.298, linetype="dashed")+ annotate("text", label = "Mean Distance (2598 bp)", size = 2, x = 6000, y = 150000)

plot_b <- inter_snp_dist  %>% filter(snp_dist_bin <=100) %>% ggplot(aes(snp_dist_bin)) + geom_histogram(bins=20) + theme_bw()+ labs(x="Inter-SNP Distance (bp)", y="Count")+theme(axis.title = element_text(face="bold", size=10, colour="black"), axis.text = element_text(colour="black", size=8))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



####MAF DISTRIBUTION####

mean(snp_freq_chr$MAF)
#0.1567062


maf_snp_distribution <- snp_freq_chr %>% ggplot(aes(MAF)) + geom_histogram(bins=30) + theme_bw()+ labs(x="Minor Allele Frequency (MAF)", y="Count")+theme(axis.title = element_text(face="bold", size=10, colour="black"), axis.text = element_text(colour="black", size=8))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_vline(xintercept=0.1567062, linetype="dashed")+ annotate("text", label = "Mean MAF (0.157)", size = 2, x = 0.21, y = 44000)

snp_stats <- snp_chromo / (plot_a | plot_b) / maf_snp_distribution + plot_annotation(tag_levels = 'A')

ggsave("../figures/snp_stats.pdf", plot = snp_stats)

