#Instructions: edit the corresponding paths to match your directory paths before running code; if on mac run L5 to troubleshoot PLINK2; others start on L8

##to get plink to run on the mac, needed to delete the mac quarantine by using the following code in terminal:
xattr -d com.apple.quarantine ~/path/to/PLINK_Files/plink2_mac_20260228/plink2

##data pre processing in bcftools in terminal
###setting working directory to use bcftools
cd /path/bcftools

###creating an index file for the vcf, which is needed to run the subsetting
./bcftools index /path/2_TB_Working_Files/joint_call.LAO03M.20251202.vcf.gz

###subsetting the vcf to just contain the autosomes, the original file also has the sex chromosomes, mtDNA, and unplaced fragments
bcftools % ./bcftools view -r A1_RagTag,A2_RagTag,A3_RagTag,B1_RagTag,B2_RagTag,B3_RagTag,B4_RagTag,C1_RagTag,C2_RagTag,C3_RagTag,D1_RagTag,D2_RagTag,D3_RagTag,D4_RagTag,E1_RagTag,E2_RagTag,E3_RagTag \
/path/1_Data/1_Working_Files/joint_call.LAO03M.20251202.vcf.gz \
-O z -o /path/1_Data/1_Working_Files/joint_call_autosomes.vcf.gz

###removing indels and non biallelic snps from the vcf files in bcftools
####this code first normalizes the data set, separating the individual alleles. This allows the code to then remove indels (as some indels can be an alt allele to a snp),
####then it removes any sites that have only 1 allele, and any that the have more than two -- essentially the biallelic filter
####in total, this allows the indel filter to catch all instances of indels, then remove uninformative sites and multiallelic sites leaving only biallelic snps
./bcftools norm -m -any \ /path/1_Data/1_Working_Files/Filter_testing/joint_call_autosomes.vcf.gz \ | ./bcftools view -v snps -m2 -M2 \
-O z -o /path/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz
          
###filtering for genotype quality scores in bcftools -- score of 9 used as min
####code tells bcftools to exclude (-e) any snps that have a genotype score or less than 9, then reads in the file, and outputs the new filtered file
./bcftools filter -e 'FMT/GQ < 9' \ /path/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz \ -O z -o allInd_SNPs_autosomes_bi_gq9.vcf.gz