PLINKpath<-"F:/2_TB_Working_Files/Plink_files/WindowsPLINK"
dir.create("./data/processed/forNe")
setwd("./data/processed/forNe")

####SNP filtering to create Ne data set for all wild LEPA####
#no locus missingness: --geno 0
#maf: 0.05
###Wild
system(paste0(PLINKpath, "/plink --bfile ../SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --maf 0.05 --geno 0 --hwe 1e-6 --make-bed --out wild_Ne_no_miss"))
#Total genotyping rate is 0.883688.
#2935647 variants removed due to missing genotype data (--geno).
#--hwe: 7 variants removed due to Hardy-Weinberg exact test.
#9506 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#5191 variants and 44 samples pass filters and QC.

####Export full vcf
system(paste0(PLINKpath, "/plink2 --bfile wild_Ne_no_miss --chr-set 17 --export vcf-4.2 bgz --out wild_Ne_no_miss"))

##split into ranch and refuge for population level nucleotide diversity
system(paste0(PLINKpath,"/plink --bfile wild_Ne_no_miss --keep ../../inputs/ranch_subset_vcftools.txt --chr-set 17 --make-bed --out ranch_wild_Ne_no_miss"))#ranch
system(paste0(PLINKpath, "/plink --bfile wild_Ne_no_miss --remove ../../inputs/ranch_subset_vcftools.txt --chr-set 17 --make-bed --out refuge_wild_Ne_no_miss")) #refuge

#export subset files
system(paste0(PLINKpath,"/plink2 --bfile ranch_wild_Ne_no_miss --chr-set 17 --export vcf-4.2 bgz --out ranch_wild_Ne_no_miss"))
system(paste0(PLINKpath,"/plink2 --bfile refuge_wild_Ne_no_miss --chr-set 17 --export vcf-4.2 bgz --out refuge_wild_Ne_no_miss"))

#split by source pop
##individuals that grouped with Ranch or Refuges in PCA and admixture were placed in their source population
##admixed individuals stayed with the population they are currently living in
#use input ranch_subset_vcftools_source_pop.txt
##split into ranch and refuge for population level nucleotide diversity
system(paste0(PLINKpath,"/plink --bfile wild_Ne_no_miss --keep ../../inputs/ranch_subset_vcftools_source_pop.txt --chr-set 17 --make-bed --out ranch_source_wild_Ne_no_miss"))#ranch
system(paste0(PLINKpath, "/plink --bfile wild_Ne_no_miss --remove ../../inputs/ranch_subset_vcftools_source_pop.txt --chr-set 17 --make-bed --out refuge_source_wild_Ne_no_miss")) #refuge

#export subset files
system(paste0(PLINKpath,"/plink2 --bfile ranch_source_wild_Ne_no_miss --chr-set 17 --export vcf-4.2 bgz --out ranch_source_wild_Ne_no_miss"))
system(paste0(PLINKpath,"/plink2 --bfile refuge_source_wild_Ne_no_miss --chr-set 17 --export vcf-4.2 bgz --out refuge_source_wild_Ne_no_miss"))



####MAF 0.01 only kept 5583 which is only 392 more SNPs...run with GONE2 and currentNe with whole dataset; if not different go with 0.05
#no locus missingness: --geno 0
#maf: 0.01
###Wild
system(paste0(PLINKpath, "/plink --bfile ../SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --maf 0.01 --geno 0 --hwe 1e-6 --make-bed --out wild_Ne_no_miss_maf01"))
#Total genotyping rate is 0.883688.
#2935647 variants removed due to missing genotype data (--geno).
#--hwe: 7 variants removed due to Hardy-Weinberg exact test.
#9114 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#5583 variants and 44 samples pass filters and QC.

####Export full vcf
system(paste0(PLINKpath, "/plink2 --bfile wild_Ne_no_miss_maf01 --chr-set 17 --export vcf-4.2 bgz --out wild_Ne_no_miss_maf01"))


###no MAF fitler
#no locus missingness: --geno 0
###Wild
system(paste0(PLINKpath, "/plink --bfile ../SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --geno 0 --hwe 1e-6 --make-bed --out wild_Ne_no_miss_nomaf"))
#Total genotyping rate is 0.883688.
#2935647 variants removed due to missing genotype data (--geno).
#--hwe: 7 variants removed due to Hardy-Weinberg exact test.
#14697 variants and 44 samples pass filters and QC.

####Export full vcf
system(paste0(PLINKpath, "/plink2 --bfile wild_Ne_no_miss_nomaf --chr-set 17 --export vcf-4.2 bgz --out wild_Ne_no_miss_nomaf"))



####filter wild_kin_roh_filter
#ROH and kinship filters; no MAF, miss 90, biallelic, coverage depth 7, genotype quality 9
#add --mac 1 to filter out polymorphic loci since GONE2 already autofilters them, but I think something weird is going on
system(paste0(PLINKpath,"/plink --bfile ../SNP_wild_dp7_gq9_bi --chr-set 17 --keep-allele-order --geno 0.1 --mac 1 --hwe 1e-6 --make-bed --out wild_kin_roh_filter_mac1"))
#Total genotyping rate is 0.883688.
#Total genotyping rate is 0.883688.
#1691187 variants removed due to missing genotype data (--geno).
#--hwe: 566 variants removed due to Hardy-Weinberg exact test.
#781354 variants removed due to minor allele threshold(s)
#477244 variants and 44 samples pass filters and QC.

####Export full vcf
system(paste0(PLINKpath, "/plink2 --bfile wild_kin_roh_filter_mac1 --chr-set 17 --export vcf-4.2 bgz --out wild_kin_roh_filter_mac1"))
