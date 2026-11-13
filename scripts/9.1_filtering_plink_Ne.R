PLINKpath<-"path/WindowsPLINK"
dir.create("./data/processed/forNe")
setwd("./data/processed/forNe")

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


##split into ranch and refuge for genetic population level
#this sorts them based on how they were grouped in the PCA: i.e. LO03M, OM331, E32M, E29M were all grouped as Ranch
##even thought they were captured in Refuge; admixed individuals stayed Ranch (E35M) or Refuge (LO01F) based on capture location
system(paste0(PLINKpath,"/plink --bfile wild_kin_roh_filter_mac1 --keep ../../inputs/ranch_subset_vcftools_source_pop.txt --chr-set 17 --make-bed --out wild_kin_roh_filter_mac1_ranch"))#ranch
system(paste0(PLINKpath, "/plink --bfile wild_kin_roh_filter_mac1 --remove ../../inputs/ranch_subset_vcftools_source_pop.txt --chr-set 17 --make-bed --out wild_kin_roh_filter_mac1_refuge")) #refuge

#export subset files
system(paste0(PLINKpath,"/plink2 --bfile wild_kin_roh_filter_mac1_ranch --chr-set 17 --export vcf-4.2 bgz --out wild_kin_roh_filter_mac1_ranch"))
system(paste0(PLINKpath,"/plink2 --bfile wild_kin_roh_filter_mac1_refuge --chr-set 17 --export vcf-4.2 bgz --out wild_kin_roh_filter_mac1_refuge"))
