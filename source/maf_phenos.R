#calculate MAF for each traits genotype file.

#load phenotype names
pheno_list <- read_csv("outputs/pheno_list.txt", 
                       col_names = FALSE)
pheno_list <- pheno_list$X1

#set paths to files.
input="/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/"
output="/project/6003429/myles_lab/abc_gwas/big_gwas_analysis/gwas/genotype_data/maf/"

#make raw file to calculate the effective number of markers.
geno_raw_list = list()
for (i in pheno_list){
  script=paste0("/project/6003429/myles_lab/bin/plink_linux_x86_64/plink --file ",input,i,"_geno_filtered --out ",output,i,"_geno_filtered --freq")
  geno_raw_list[[i]] = script
}
write.table(paste(geno_raw_list[1:38]), "shell_scripts/maf_phenos.sh", sep="\t", row.names = F, quote = F, col.names = F)
