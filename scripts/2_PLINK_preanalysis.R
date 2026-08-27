#Pre-Analysis formatting with PLINK
##Common pitfalls: 1) wrong PLINK binary called (i.e. trying to run MAC-OS executable on Windows/Linux, 
##2) check your paths, 3) check how your paths are printed if cloud storage added pesky spaces and & symbols)

#Define paths
##update with your path to PLINK and PLINK2
PLINKpath<-"F:/2_TB_Working_Files/Plink_files/WindowsPLINK"

##inputs on Zenodo (ADD DOI): allInd_SNPs_autosomes.vcf.gz
##update with your path
ZenodoVCFpath<-"F:/2_TB_Working_Files/Manu_Files/allInd_SNPs_autosomes.vcf.gz"

##output path
OUT<- paste0(getwd(), "/data/processed/")
OUT_name<-shQuote(paste0(OUT, "SNP_AllChrom_AllInd_dp7_gq9_bi"), type = "cmd") #change type as appropriate for your OS shell

#PLINK pre-analysis setup
##creating plink files from the base file created from BCFtools
system(paste0(PLINKpath,"/plink2 --vcf ", ZenodoVCFpath," --keep-allele-order --allow-extra-chr --vcf-min-dp 7 --vcf-max-dp 22 --max-alleles 2 --chr-set 17 --make-bed --out ", OUT_name))
###89 individuals 103792041 variants remain after filter for depth and biallelic

##creating base plink files pre-standard filtering
###renaming chromosomes to a standard number based, which is what plink is expecting
bim <- read.table("./data/processed/SNP_AllChrom_AllInd_dp7_gq9_bi.bim", stringsAsFactors = FALSE)
colnames(bim) <- c("chr", "snp", "cm", "pos", "a1", "a2")
###get the unique chr names
unique_chrs <- unique(bim$chr)
print(unique_chrs)
chr_map <- data.frame(
  old_chr = unique_chrs,
  new_chr = 1:length(unique_chrs),
  stringsAsFactors = FALSE
)
print(chr_map)
###apply that to the bim file directly
for (i in 1:nrow(chr_map)) {
  bim$chr[bim$chr == chr_map$old_chr[i]] <- chr_map$new_chr[i]
}
###write updated bim file
write.table(bim, "SNP_AllInd_unfilt_chrfix.bim", quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)

##have plink write new bim bam bed files with the new bim file we created
system("./plink --bed SNP_AllChrom_AllInd_dp7_gq9_bi.bed --bim SNP_AllInd_unfilt_chrfix.bim --fam SNP_AllChrom_AllInd_dp7_gq9_bi.fam --make-bed --out SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix")

###give snp ids -- labels every snps with a unique code
system("./plink2 --bfile SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix --chr-set 17 --set-all-var-ids @_# --new-id-max-allele-len 20 --make-bed --out SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix_uniqueID")

###subset data -- remove mountain lion and duplicates, for manuscript 1 keeping only wild populations
system("./plink --bfile SNP_AllChrom_AllInd_dp7_gq9_bi_chrfix_uniqueID --keep wild_subset.txt --chr-set 17 --make-bed --out SNP_wild_dp7_gq9_bi")
#Total genotyping rate in remaining samples is 0.861001; 103792041 variants and 44 samples pass filters and QC.
#at this point, coverage depth of 7, genotype quality of 9, and biallelic filters applied to create base file

####SNP filtering to create standard data set for all wild LEPA####
#apply filters -- maf, miss, hwe> base, other adjustments can follow
###Wild
system("./plink --bfile SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --maf 0.05 --geno 0.1 --hwe 1e-6 --make-bed --out wild_standard_final")
#66269934 variants removed due to missing genotype data (--geno)
#3576 variants removed due to Hardy-Weinberg exact test
#33435588 variants removed due to minor allele threshold(s)
#4082943 variants and 44 samples pass filters and QC.

####Export full wild standard vcf
system("./plink2 --bfile wild_standard_final --chr-set 17 --export vcf-4.2 bgz --out wild_standard_final")
##split into ranch and refuge for population level nucleotide diveristy
system("./plink --bfile wild_standard_final --keep ranch_subset_vcftools.txt --chr-set 17 --make-bed --out ranch_standard_postsubset") #ranch
system("./plink --bfile wild_standard_final --remove ranch_subset_vcftools.txt --chr-set 17 --make-bed --out refuge_standard_postsubset") #refuge
#export subsetted files
system("./plink2 --bfile ranch_standard_postsubset --chr-set 17 --export vcf-4.2 bgz --out ranch_standard_postsubset")
system("./plink2 --bfile refuge_standard_postsubset --chr-set 17 --export vcf-4.2 bgz --out refuge_standard_postsubset")

####additional filtering for select analyses####
#LD pruning
system("./plink --bfile wild_standard_final --chr-set 17 --keep-allele-order --indep-pairwise 50 5 0.5 --out wild_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
#change to indep-pairwise instead of indep
system("./plink --bfile wild_standard_final --extract wild_LDpruned_0.5_out.prune.in --chr-set 17 --make-bed --out wild_LDpruned_05") #extract SNPs and create new files
#Total genotyping rate is 0.925318; 273852 variants and 44 samples pass filters and QC..
#write vcf
system("./plink2 --bfile wild_LDpruned_05 --chr-set 17 --export vcf-4.2 bgz --out wild_LDpruned_05")


#ROH and kinship filters; no MAF, miss 90, biallelic, coverage depth 7, genotype quality 9
#also used in kinship analyses; as is recommended not to filter for MAF or LD prune
system("./plink --bfile SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --geno 0.1 --hwe 1e-6 --make-bed --out wild_kin_roh_filter")
#for bcftools roh selection: ld pruning is required as it assumes every base is independent
system("./plink --bfile wild_kin_roh_filter --chr-set 17 --keep-allele-order --indep-pairwise 50 5 0.5 --out roh_LDpruned_0.5_out") #makes an out and in files of SNps to keep and SNPs to remove
#change to indep-pairwise instead of indep
system("./plink --bfile wild_kin_roh_filter --extract roh_LDpruned_0.5_out.prune.in --chr-set 17 --make-bed --out roh_LDpruned_05") #extract SNPs and create new files
#Total genotyping rate is 0.925318; 273852 variants and 44 samples pass filters and QC..
#write vcf
system("./plink2 --bfile roh_LDpruned_05 --chr-set 17 --export vcf-4.2 bgz --out roh_LDpruned_05")
#2216544 variants and 44 samples pass filters and QC.
#output as vcf for use in bcftools
system("./plink2 --bfile wild_kin_roh_filter --chr-set 17 --export vcf-4.2 bgz --out wild_kin_roh_filter")
