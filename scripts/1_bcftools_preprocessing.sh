#Orginal author: Tyler Bostwick
#Edited: Ariana Cerreta, 17-Sept-2026

#Instructions: edit the corresponding paths to match your directory paths before running code; if on mac run L10 to troubleshoot PLINK2; others start on L13
##if you are opening this code in R from a forked repository, stay in default project directory from tb_ocelot_genetics.Rproj
##you will need to download joint_call.LAO03M.20251202.vcf.gz from Zenodo and know what directory it is in
##you will need to download ocelot genome of ocelot genome GCA_058742865.1 as a .fna.gz from NCBI
##TIP: to run these lines within an R session, click on the "Terminal" in the console pane

##to get plink to run on the mac, needed to delete the mac quarantine by using the following code in terminal:
xattr -d com.apple.quarantine ~/path/to/PLINK_Files/plink2

##data pre processing in bcftools in terminal
#paths (edit according your directory set up)
export BCFTOOLS="/path/bcftools-1.24/bin" #directory to your executable bcftools 
export GENOME="/path/GCA_058742865.1_TAMU_Opar_alb_1.0_genomic.fna" #directory ending in .fna.gz of ocelot genome GCA_058742865.1
export CHROM_MAP="./data/inputs/TAMU_NCBI_Chrom_key_LEPA_ocelot.txt" #you should not need to edit if forked properly from GitHub
export INPUT_FOLDER="/path" #directory location where you placed joint_call.LAO03M.20251202.vcf.gz from Zenodo
export OUTPUT_FOLDER="/path/file/write/directory" #where all the files created from original file are written; you may have to make this

###creating an index file for the vcf, which is needed to run the subsetting
"${BCFTOOLS}"/bcftools index "${INPUT_FOLDER}"/joint_call.LAO03M.20251202.vcf.gz

###subsetting the vcf to just contain the autosomes, the original file also has the sex chromosomes, mtDNA, and unplaced fragments
"${BCFTOOLS}"/bcftools view -r A1_RagTag,A2_RagTag,A3_RagTag,B1_RagTag,B2_RagTag,B3_RagTag,B4_RagTag,C1_RagTag,C2_RagTag,C3_RagTag,D1_RagTag,D2_RagTag,D3_RagTag,D4_RagTag,E1_RagTag,E2_RagTag,E3_RagTag \
"${INPUT_FOLDER}"/joint_call.LAO03M.20251202.vcf.gz \
-O z -o "${OUTPUT_FOLDER}"/joint_call_autosomes.vcf.gz

###rename the chromosomes in the vcf so that the chromosome names match the public repository genome
"${BCFTOOLS}"/bcftools annotate --rename-chrs "$CHROM_MAP" "${OUTPUT_FOLDER}"/joint_call_autosomes.vcf.gz \
-O z -o "${OUTPUT_FOLDER}"/joint_call_autosomes_NCBIchrom.vcf.gz

###removing indels and non biallelic snps from the vcf files in bcftools
####this code first normalizes the data set using the reference genome (.fna.gz), removes duplicates, and keeps all variant types merged
####then it only keeps SNPs (1 bp variant) removes any sites that have only 1 allele, and any that the have more than two -- essentially the biallelic filter
####successfully removes multiallelic, multitype loci like CM178980.1       29995796        .       GGTACCCAAGT     AGTACCCAAGT,G
####successfully removes SNP-type MNPS like CM178977.1       1236443 .       G       T,A

"${BCFTOOLS}"/bcftools norm --fasta-ref "$GENOME" -d exact --multiallelics +any \
"${OUTPUT_FOLDER}"/joint_call_autosomes_NCBIchrom.vcf.gz \ | "${BCFTOOLS}"/bcftools view -v snps -m2 -M2 \
-O z -o "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes.vcf.gz

####Optional Check####
"${BCFTOOLS}"/bcftools index "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes.vcf.gz
"${BCFTOOLS}"/bcftools index -n "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes.vcf.gz #should be 94,736,326
"${BCFTOOLS}"/bcftools view -H -r CM178977.1:1236443 "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes.vcf.gz #should have no output
"${BCFTOOLS}"/bcftools view -H -r CM178980.1:29995796  "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes.vcf.gz #should have no output
####End Optional Check####

###filtering for genotype quality scores in bcftools -- score of 9 used as min
####code tells bcftools to exclude (-e) any snps that have a genotype score or less than 9, then reads in the file, and outputs the new filtered file
"${BCFTOOLS}"/bcftools filter -e 'FMT/GQ < 9' "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes.vcf.gz \
-O z -o "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes_bi_gq9.vcf.gz

####Optional Check####
"${BCFTOOLS}"/bcftools index "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes_bi_gq9.vcf.gz
"${BCFTOOLS}"/bcftools index -n "${OUTPUT_FOLDER}"/allInd_SNPs_autosomes_bi_gq9.vcf.gz #should be 2,950,351
