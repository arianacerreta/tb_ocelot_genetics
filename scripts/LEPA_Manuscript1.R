#04/22/26 - T. Bostwick
#data processing and analyses for Manuscript 1: wild ocelot genomics
#data aligned to the new ocelot genome (Foley et al. 2026; https://doi.org/10.1093/jhered/esaf095)
#vcf from the lab had minimal filtering, only basic quality control, and no imputation

###~~~ Filtering steps~~~
#First in BCFTools:
    #select only autosomes
    #remove indels and non biallelic snps
    #filter for a genotype quality score above 9
#next, in plink for processesing ease:
    #rename chromosomes
    #give snps unique id's
    #for manuscript 1 -- subset just the wild
#then, filtering in plink for each analysis
    #depth of 7
    #generally, maf, miss, and hwe filters applied here

####setting up####
#setting working directories
setwd("/Volumes/Expansion/2_TB_Working_Files/Plink_files") #WD if working out of the hard drive
setwd("~/Documents/Masters_Work/Analyses/1_Data/1_Working_Files") #WD if working off of the laptop
#installing and libraring required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(ggplot2)
library(dplyr)
library(HardyWeinberg)
library(snpStats)
library(adegenet)
library(poppr)
library(vcfR)
library(promises)
library(adegenet)
library(vcfR)
library(poppr)
library(SNPRelate)
library(dartR)
##for karyotype plots
install.packages(c("karyoploteR", "regioneR", "GenomicRanges", "data.table", "IRanges", "GenomicAlignments"))
BiocManager::install("karyoploteR")
if (!requireNamespace("regioneR", quietly = TRUE))
  BiocManager::install("regioneR")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")
if (!requireNamespace("IRanges", quietly = TRUE))
  BiocManager::install("IRanges")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
library(karyoploteR)
library(regioneR)
library(GenomicRanges)
library(data.table)
library(IRanges)
library(GenomicAlignments)
##for admixture plots
library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)
library(tidyr)
library(stringr)
library(dplyr)

#troubleshooting code for mac
#to get plink to run on the mac, needed to delete the mac quarantine by using the following code in terminal:
    #xattr -d com.apple.quarantine ~/Documents/Masters_Work/Analyses/PLINK_Files/plink2_mac_20260228/plink2


################################################################################
####data pre processing in bcftools in terminal####
#setting working directory to use bcftools
  #cd /Users/tylerbostwick/bcftools
#creating an index file for the vcf, which is needed to run the subsetting
  #./bcftools index /Volumes/Expansion/2_TB_Working_Files/joint_call.LAO03M.20251202.vcf.gz
#subsetting the vcf to just contain the autosomes, the original file also has the sex chromosomes, mtDNA, and unplaced fragments
  #bcftools % ./bcftools view -r A1_RagTag,A2_RagTag,A3_RagTag,B1_RagTag,B2_RagTag,B3_RagTag,B4_RagTag,C1_RagTag,C2_RagTag,C3_RagTag,D1_RagTag,D2_RagTag,D3_RagTag,D4_RagTag,E1_RagTag,E2_RagTag,E3_RagTag \
  #/Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/joint_call.LAO03M.20251202.vcf.gz \
  #-O z -o /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/joint_call_autosomes.vcf.gz

#removing indels and non biallelic snps from the vcf files in bcftools
  #./bcftools norm -m -any \ /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/joint_call_autosomes.vcf.gz \ | ./bcftools view -v snps -m2 -M2 \
  #-O z -o /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz
          #this code first normalizes the data set, separating the individual alleles. This allows the code to then remove indels (as some indels can be an alt allele to a snp),
          #then it removes any sites that have only 1 allele, and any that the have more than two -- essentially the biallelic filter
          #in total, this allows the indel filter to catch all instances of indels, then remove uninformative sites and multiallelic sites leaving only biallelic snps

#filtering for genotype quality scores in bcftools -- score of 9 used as min
#./bcftools filter -e 'FMT/GQ < 9' \ /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz \ -O z -o allInd_SNPs_autosomes_bi_gq9.vcf.gz
      #code tells bcftools to exclude (-e) and snps that have a genotype score or less than 9, then reads in the file, and outputs the new filtered file

####creating plink files from the base file created from BCFtools####
system("./plink2 --vcf allInd_SNPs_autosomes.vcf.gz --keep-allele-order --allow-extra-chr --vcf-min-dp 7 --vcf-max-dp 22 --max-alleles 2 --chr-set 17 --make-bed --out SNP_AllChrom_AllInd_dp7_gq9_bi")
      #89 individuals 103792041 variants remain after filter for depth and biallelic

####creating base plink files pre-standard filtering####
##renaming chromosomes to a standard number based, which is what plink is expecting
bim <- read.table("SNP_AllChrom_AllInd_dp7_gq9_bi.bim", stringsAsFactors = FALSE)
colnames(bim) <- c("chr", "snp", "cm", "pos", "a1", "a2")
#get the unique chr names
unique_chrs <- unique(bim$chr)
print(unique_chrs)
chr_map <- data.frame(
  old_chr = unique_chrs,
  new_chr = 1:length(unique_chrs),
  stringsAsFactors = FALSE
)
print(chr_map)
#apply that to the bim file directly
for (i in 1:nrow(chr_map)) {
  bim$chr[bim$chr == chr_map$old_chr[i]] <- chr_map$new_chr[i]
}
#write updated bim file
write.table(bim, "SNP_AllInd_unfilt_chrfix.bim", quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)
#have plink write new bim bam bed files with the new bim file we created
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

##############################################################
####Data Analysis
####PCA Analysis -- Done!####
##PCA on Wild individuals --- done!
#reading in data
system("./plink --bfile wild_standard_final --pca --chr-set 17 --out pca_wild") #run the pca code, specified the number of chromosomes
pca.data.wild <- read.table("pca_wild.eigenvec", header=FALSE) #load pca results from the eigenvec file
colnames(pca.data.wild)[2] <- "ID" #rename second column to allow for appending pop assignment later
eigenvalues <- read.table("pca_wild.eigenval", header=FALSE)$V1 #read in eigenvalues to calc variance
#calc variance
total_variance <- sum(eigenvalues)
pc1_variance <- round((eigenvalues[1] / total_variance) * 100, 2)
pc2_variance <- round((eigenvalues[2] / total_variance) * 100, 2)
#adding population information
wild_origins <- read.csv("wild_origins.csv", header = TRUE) #read in origins file
pca.data.wild.origins <- left_join(pca.data.wild, wild_origins, by = "ID") #append origin data to pca data
pca.data.wild.origins <- pca.data.wild.origins %>%
  mutate(ID = gsub("-.*", "", ID))
#plotting
wild_pca <- ggplot(pca.data.wild.origins, aes(x=V3,y=V4)) +  #plot with individual ID's and by origin
  geom_point(aes(shape = Pop, color = Pop), size = 3) +
  geom_rect(data = subset(pca.data.wild.origins, ID %in% c("E35M", "LO01F")),
            aes(xmin = min(V3) - 0.02, xmax = max(V3) + 0.04,
                ymin = min(V4) - 0.02, ymax = max(V4) + 0.04),
            fill = NA, color = "black", linewidth = 1, linetype = "dashed") +
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("LO03M")), #in correct position
            aes(label=ID), vjust=-0.6, hjust=-0.1, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("OM331")), #in correct position
            aes(label=ID), vjust=-0.6, hjust=-0.1, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("E35M")), #in correct position
            aes(label=ID), vjust=-0.6, hjust=-0.1, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("LO01F")), #in correct position
            aes(label=ID), vjust=-1, hjust=0.02, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("E29M")), #in correct position
            aes(label=ID), vjust=1.5, hjust=-0.1, size=4, color = "black") + 
  geom_text(data = subset(pca.data.wild.origins, ID %in% c("E32M")),  #in correct position
            aes(label=ID), vjust=-0.4, hjust=-0.2, size=4, color = "black") + 
  scale_color_manual(name = "Origin", values = c("Refuge" = "#01004c", "Ranch" = "orchid")) +
  scale_shape_manual(name = "Origin", values = c("Refuge" = 19, "Ranch" = 17)) +
  labs(x = paste0("PC1 (", pc1_variance, "%)"),
       y = paste0("PC2 (", pc2_variance, "%)")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
wild_pca
ggsave("wild_pca.png", dpi = 300)

####Kinship -- Done!####
##king-robust kingship estimator for WILD individuals 
system("./plink2 --bfile wild_kin_roh_filter --make-king-table --allow-extra-chr --chr-set 17 --out wild_KING_manu")
king.wild.matrix <- read.table("wild_KING_manu.kin0", header = FALSE) #make king table into object
colnames(king.wild.matrix) <- c("#FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header
#remove "-x" from the individual id's
king.wild.matrix <- king.wild.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
#make csv from the .kin0 file
write.csv(king.wild.matrix, "wild_pairwise_kinship_manu.csv")

#heat map of kinship for wild individuals -- binned
ggplot(data = king.wild.matrix, aes(x=IID1, y=IID2, fill = KINSHIP)) +
  geom_tile(color="white") +
  scale_fill_stepsn(name = "Kinship", breaks = c(0, 0.04, 0.1, 0.2),
                    limit = c(0, 0.3),
                    labels = c("Unrelated", "3rd Degree", "2nd Degree", "1st Degree"),
                    aesthetics = "fill",
                    space = "Lab",
                    colors = c("grey100", "gold1", "orangered", "firebrick"))+
  labs(x="Individuals", y="Individuals")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

####Admixture -- Done!####
##admixture plotting
#using LD pruned dataset
#.q files from ADMIXTURE performed on a mac computer; code as follows:

#setting working directory to folder with admixture
  #cd /Users/tylerbostwick/Documents/Masters_Work/Analyses/ADMIXTURE/dist/admixture_macosx-1.3.0
#running admixture
  #for k in {2..5}; do ./admixture --cv=10 wild_LDpruned_05.bed $k > wild_LDpruned_05k${k}.txt; done
    #runs the admixture program for k 2-5, generating cross-validation, and writing it to a txt file for each k ran

###wild admixture plot
#reading in and preparing data
w_k2_table <- read.table("wild_LDpruned_05.2.q") ## finshed through here
w_ids <- read.table("wild_pop_samp_id.txt")
w_ids = rename(w_ids, sample_id = V1, pop_id = V2)
w_k2 <- cbind(w_ids, w_k2_table)
w_k2_long <- pivot_longer(
  w_k2, cols = c(V1, V2),
  names_to = "ancestry",
  values_to = "proportion"
) #pivots the table to the proportions are able to be plotted as stacked bars
w_k2_long$pop_id <- str_to_upper(w_k2_long$pop_id)
#plot code
w_k2plot <-
  ggplot(w_k2_long, aes(x = factor(sample_id), y = proportion, fill = ancestry)) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(pop_id), switch = "both", scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, size = 15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 12),
        strip.text.x = element_text(face = "bold", size = 15),
        strip.placement = "outside",
        axis.text.y = element_text(size = 15),
  ) +
  scale_fill_manual(values = c("V1" = "#01004c", "V2" = "#ffb2b0")) +
  guides(fill = "none")
w_k2plot
ggsave("wild_k2_admixture_Apr2026.png", w_k2plot, width = 15, height = 8, bg = "white")
####ROH -- Done!####
#set working directory to roh specific folder
setwd("~/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/manu_roh_selection")


##experirmental -- used bcftools to select roh, uses a HMM to select runs
#code run in termainl is as follows:

#clean the header of the vcf for use in bcftools- removes a flag placed by plink that causes errors
#./bcftools view -h /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/manu_roh_selection/roh_LDpruned_05.vcf.gz | grep -v "##chrSet" > clean_header.txt
#./bcftools reheader -h clean_header.txt /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/manu_roh_selection/roh_LDpruned_05.vcf.gz -o /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/manu_roh_selection/clean_roh_LDpruned_05.vcf.gz

#identify roh:
#./bcftools roh -G30 --estimate-AF - --rec-rate 1.1e-8 -Or -o roh_LD05.txt \ /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/manu_roh_selection/clean_roh_LDpruned_05.vcf.gz
    #key changes, now estimates allele frequencies from data instead of using default, uses the domestic cat recombination rate
    #and is using LD pruned data to better fit model assumptions
    #as of 4/29/26 -- fixed the rec-rate flag
        #kept g30, lowered it down to g20 and had minimal change in percentages, no break in long runs with the correct rec rate
        

###brief view of the output and froh by individual:
#read in output table, only the RG (roh segment) lines
roh_LD <- read.table("roh_LD05_v2.txt", comment.char = "#", header = FALSE)
# Keep only ROH segment rows (RG), drop per-site rows (ST)
roh_LD <- roh_LD[roh_LD$V1 == "RG", ]
# Name the columns
colnames(roh_LD) <- c("type", "sample", "chromosome", "start", "end", "length_bp", "n_markers", "quality")
# Drop the type column since everything is RG now
roh_LD$type <- NULL
##filter roh segments for high quality scores only
# Keep only high confidence ROH
roh_LD_hq <- roh_LD[roh_LD$quality >= 30, ]
# Filter by minimum length (e.g. 500kb, similar to the plink parameters)
roh_LD_filt <- roh_LD_hq[roh_LD_hq$length_bp >= 500000, ]
##testing output: FROH by individual
#Sum ROH length per individual
ind_roh_LD <- aggregate(length_bp ~ sample, data = roh_LD_filt, FUN = sum)
# Calculate FROH
ind_roh_LD$FROH <- (ind_roh_LD$length_bp / 2314601649) * 100 #updated number of bases to reflect the geofferys cat genome

#writing tables
write.csv(ind_roh_LD, "froh_by_individual_hmm.csv")
write.csv(roh_LD_filt, "roh_segments_filt_qual_length_hmm.csv", row.names = FALSE)

#adding population labels to data
pop_id <- read.table("wild_pop_id.txt", header = FALSE)
colnames(pop_id) <- c("sample", "pop")
as.data.frame(pop_id)
wild_roh_pop <- left_join(ind_roh_LD, pop_id, by = "sample")

####ROH assessment
#average % genome in roh by population
population_mean_roh <- wild_roh_pop %>%
  group_by(pop) %>%
  summarise(Average = mean(FROH, na.rm = TRUE),
            Count = n(),
            StdDev = sd(FROH, na.rm = TRUE),
            StdError = (sd(FROH, na.rm = TRUE)/sqrt(n())))

write.csv(population_mean_roh, "population_mean_froh_table.csv")

#violin plot of wild FROH
ggplot() +
  # violin plot
  geom_violin(data = wild_roh_pop, aes(x = pop, y = FROH, fill = pop), 
              alpha = 0.7) +
  # Add individual points
  geom_jitter(data = wild_roh_pop, aes(x = pop, y = FROH), 
              width = 0.1, alpha = 0.6, size = 3) +
  # Add population means with error bars
  geom_point(data = population_mean_roh, aes(x = pop, y = Average), 
             color = "red", size = 4, shape = 18) +
  geom_errorbar(data = population_mean_roh, 
                aes(x = pop, y = Average, 
                    ymin = Average - 2*StdDev, ymax = Average + 2*StdDev),
                color = "red", width = 0.2, size = 1) +
  scale_fill_manual(values = c("ranch" = "#ffb2b0", "refuge" = "#01004c")) +
  # Labels and theme
  labs(x = "Population", y = expression(F[ROH] ("%"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")


####proportion of ROH lengths
#read in .hom files
roh_seg_df <- read.csv("roh_segments_filt_qual_length_hmm.csv", header = TRUE)
#adding population information
pop_id <- read.table("wild_pop_id.txt", header = FALSE)
colnames(pop_id) <- c("sample", "pop")
as.data.frame(pop_id)
roh_seg_df <- merge(roh_seg_df, pop_id, by = "sample", all = TRUE)
#calculate length in Mb -- bp to mb conversion
roh_seg_df$length_MB <- roh_seg_df$length_bp/1000000
roh_seg_df$length_KB <- roh_seg_df$length_bp/1000 #converting bp to kilobase pair
#look at resulting distribution
hist(roh_seg_df$length_MB, main="Distribution of ROH lengths", xlab="Length (MB)")
max(roh_seg_df$length_MB, na.rm = TRUE) #105.4423
min(roh_seg_df$length_MB, na.rm = TRUE) #0.500084
mean(roh_seg_df$length_MB, na.rm = TRUE) #4.747192
#define length categories
roh_seg_df$Category <- cut(roh_seg_df$length_MB,
                           breaks = c(0, 1, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, Inf),
                           labels = c("<1Mb", "1-2Mb", "2-4Mb", "4-6Mb", "6-8Mb", "8-10MB", "10-20Mb",
                                      "20Mb-30Mb", "30Mb-40Mb", "40Mb-50Mb", "50Mb-60Mb", "60Mb-70Mb",
                                      "70Mb-80Mb", "80Mb-90Mb", "90Mb-100Mb", ">100Mb"),
                           include.lowest = TRUE)
#write new table
write.csv(roh_seg_df, "hmm_roh_seg_categorized.csv")
#get total proportion -- cumulative
total_length_all <- sum(roh_seg_df$length_MB)
total_summary <- roh_seg_df %>%
  group_by(Category) %>%
  summarize(
    Count = n(),
    Total_Length_MB = sum(length_MB),
    Proportion_Length = sum(length_MB)/total_length_all,
    Proportion_count = n()/nrow(roh_seg_df)
  )
#save results
write.csv(total_summary, "overall_roh_length_proportions.csv")

#summary by individual
ind_roh_summary <- roh_seg_df %>% 
  group_by(sample) %>% 
  summarize(
    Count = n(),
    total_length = sum(length_bp),
  )
ind_roh_summary$FROH <- (ind_roh_summary$total_length / 2314601649) * 100
sum(ind_roh_summary$Count)
mean(roh_seg_df$length_KB)

ind_bin_summary <- roh_seg_df %>% 
  group_by(sample, Category) %>% 
  summarize(
    count = n()
  )
view(ind_bin_summary)
#plot generation time to roh segment
roh_seg_cat <- read.csv("hmm_roh_seg_categorized.csv", header = TRUE)
roh_seg_cat$gen <- 100/(2*(1.1)*roh_seg_cat$length_MB) #calculating generation time

#binning generation time
roh_seg_cat <- roh_seg_cat %>% 
  mutate(
    gen_bins = case_when(
      gen < 2 ~ "<2",
      gen >= 2 & gen < 5 ~ "5-2",
      gen >= 5 & gen < 10 ~ "10-5",
      gen >= 10 & gen < 20 ~ "20-10",
      gen >= 20 & gen < 40 ~ "40-20",
      gen >= 40 & gen < 60 ~ "60-40",
      gen >= 60 & gen < 80 ~ "80-60",
      gen >= 12 ~ "80+",
      TRUE ~"Other"
    ),
    gen_bins = factor(gen_bins, levels = c("80+", "80-60", "60-40", "40-20", "20-10", "10-5", "5-2", "<2"))
  )
#get ranges of roh length in each generaton bin
bin_ranges <- roh_seg_cat %>% 
  group_by(gen_bins) %>% 
  summarise(
    min_length = min(length_MB, na.rm =TRUE),
    max_length = max(length_MB, na.rm = TRUE),
    range_span = max_length - min_length,
    .groups = 'drop'
  )
#merge bin ranges with data frame
gen_plot_df <- roh_seg_cat %>% 
  left_join(bin_ranges, by = "gen_bins") %>% 
  mutate(range_label = paste0(min_length, "-", max_length, " MB"))
#plotting generation time by population
ggplot(gen_plot_df, aes(x = gen_bins, fill = pop)) +
  geom_bar(position = "dodge", stat = "count") +
  scale_y_continuous(breaks = breaks_width(100)) +
  labs(
    x = "Generation Time",
    y = "ROH Count",
    fill = "Populations"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("ranch" = "orchid", "refuge" = "#01004c"),
                    labels = c("ranch" = "Ranch", "refuge" = "Refuge")) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )


#summary by population 
population_category_counts <- roh_seg_df %>%
  group_by(pop, Category)  %>%
  summarise(Count = n(), Total_Length_MB = sum(length_MB), .groups = "drop")
#proportion by population
population_category_proportions <- population_category_counts %>%
  group_by(pop) %>%
  mutate(Proportion_Count = Count/sum(Count),
         Proportion_Length = Total_Length_MB/sum(Total_Length_MB))

#normalizing data for direct comparison
individuals_per_pop <- data.frame(
  pop = c("ranch", "refuge"),
  n_individuals = c(23, 21)
)

normalized_roh <- population_category_counts %>%
  left_join(individuals_per_pop, by = "pop") %>%
  group_by(pop) %>%
  mutate(
    # Within-population proportions
    Proportion_Count = Count / sum(Count),
    Proportion_Length = Total_Length_MB / sum(Total_Length_MB),
    
    # Per-individual normalized metrics
    Avg_ROH_Count_Per_Individual = Count / n_individuals,
    Avg_ROH_Length_MB_Per_Individual = Total_Length_MB / n_individuals
  ) %>%
  ungroup()
write.csv(normalized_roh, "pop_normalized_roh.csv")


####Visualizing ROH --karyotype plot --- individual plots working, need to find and update chromosome lengths -- email brian davis?
#read in hom file
roh_seg_df <- read.csv("hmm_roh_seg_categorized.csv", header = TRUE)
#make data frame for plotting
plot_df <- data.frame(
  chr = paste0(roh_seg_df$chromosome),  # Add 'chr' prefix if needed
  start = roh_seg_df$start, #starting location of roh
  end = roh_seg_df$end, #ending location of roh
  kb = roh_seg_df$length_KB,   # ROH length in KB
  nsnp = roh_seg_df$n_markers,   # Number of SNPs in ROH
  sample_id = roh_seg_df$sample, # Individual ID
  pop = roh_seg_df$pop #population information
)
#applying chromosome names for plots
chr_map <- c(
  "1" = "chr1",
  "2" = "chr2",
  "3" = "chr3",
  "4" = "chr4",
  "5" = "chr5",
  "6" = "chr6",
  "7" = "chr7",
  "8" = "chr8",
  "9" = "chr9",
  "10" = "chr10",
  "11" = "chr11",
  "12" = "chr12",
  "13" = "chr13",
  "14" = "chr14",
  "15" = "chr15",
  "16" = "chr16",
  "17" = "chr17"
)
plot_df$chr <- chr_map[plot_df$chr]
#create custom feline genotype for karyoploteR, make a custom plot type for the 17 chr
LP_chr_sizes <- data.frame(
  chr = c(paste0("chr", 1:17)), 
  start = rep(1, 17),
  end = c(241490831,  # chr1
          170782778,  # chr2
          140550405,  # chr3
          207231548,  # chr4
          154917953,  # chr5
          149666277,  # chr6
          142417080,  # chr7
          222385042,  # chr8
          159437941,  # chr9
          154776197,  # chr10
          116981451,   # chr11
          90819563,   # chr12
          96629807,   # chr13
          95534813,   # chr14
          64119235,   # chr15
          62492039,   # chr16
          44368689   # chr17
  )
)  
feline_genome_gr <- GRanges(seqnames = LP_chr_sizes$chr, ranges = IRanges(start = LP_chr_sizes$start, end = LP_chr_sizes$end)) #make into grange object
ocel_genome <- feline_genome_gr #making the custom plot type
seqlevels(ocel_genome) <- LP_chr_sizes$chr
seqlengths(ocel_genome) <- LP_chr_sizes$end

# roh convert to GRanges object
roh_gr <- GRanges(
  seqnames = plot_df$chr,
  ranges = IRanges(start = plot_df$start, end = plot_df$end),
  kb = plot_df$kb,
  nsnp = plot_df$nsnp,
  sample_id = plot_df$sample_id) #making hom file into grange

pop_gr <- GRanges( #### skipped for now
  seqnames = plot_df$chr,
  ranges = IRanges(start = plot_df$start, end = plot_df$end),
  kb = plot_df$kb,
  nsnp = plot_df$nsnp,
  sample_id = plot_df$sample_id,
  population_id = plot_df$pop_id)

#####creating function for plotting -- individual
plot_individual_roh <- function(sample_id, output_file = NULL) {
  # Filter ROH data for specific individual
  individual_roh <- roh_gr[mcols(roh_gr)$sample_id %in% sample_id]
  # Determine if output should go to a file
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 7)
  }
  #create the plot using the custom genome
  ind_kp <- plotKaryotype(genome = ocel_genome, plot.type = 1, main = paste("ROH for", sample_id))
  #kpAddChromosomeNames(w_kp, srt = 45, cex = 0.8) #not using this line for now
  #plot individuals
  kpRect(ind_kp, 
         chr = as.character(seqnames(individual_roh)), 
         x0 = start(individual_roh), 
         x1 = end(individual_roh),
         y0 = 0, 
         y1 = 1, 
         col = "#FF000080",  # Semi-transparent red
         border = "#FF0000",
         r0 = 0, r1 = 1,
         data.panel = "ideogram")
  # Close file if opened
  if (!is.null(output_file)) {
    dev.off()
  }
  # Return the filtered data
  return(individual_roh)
}

####function for plotting -- population wide   ~~~~###skipped for now - not updated
plot_population_roh_smoothed <- function(population_id, output_file = NULL, window_size = 1e5, smooth = TRUE) {
  # Subset to individuals from the desired population
  pop_roh <- pop_gr[mcols(pop_gr)$population_id == population_id]
  
  if (length(pop_roh) == 0) {
    message("No ROH data found for population: ", population_id)
    return(NULL)
  }
  
  # Optional output to file
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 7)
  }
  
  # Plot karyotype
  pop_kp <- plotKaryotype(genome = ocel_genome, plot.type = 2, main = paste("ROH Density for", population_id))
  
  # Calculate raw coverage (how many ROHs overlap each base)
  roh_cov <- coverage(pop_roh)
  
  # Smooth and plot each chromosome
  for (chr in names(roh_cov)) {
    cov_vector <- as.numeric(roh_cov[[chr]])
    
    if (smooth) {
      # Apply rolling mean smoothing (simple moving average)
      kernel <- rep(1/window_size, window_size)
      smoothed <- stats::filter(cov_vector, kernel, sides = 2, circular = FALSE)
    } else {
      smoothed <- cov_vector
    }
    
    # Build GRanges object with smoothed values
    pos <- IRanges(start = seq_along(smoothed), width = 1)
    smoothed_gr <- GRanges(seqnames = chr, ranges = pos, score = smoothed)
    
    # Reduce resolution for plotting (optional)
    smoothed_gr <- smoothed_gr[start(smoothed_gr) %% 1000 == 0]  # sample every 1kb
    
    # Plot line track
    kpLines(pop_kp, chr = as.character(seqnames(smoothed_gr)), x = start(smoothed_gr), y = mcols(smoothed_gr)$score, col = "red", r0 = 0.5, r1 = 0.8)
  }
  
  if (!is.null(output_file)) {
    dev.off()
  }
  
  invisible(NULL)
}

##using the function -- plotting individual roh
unique_samples <- unique(plot_df$sample_id)

##using the function -- plotting population roh. ~~~~~##skipped for now, not updated
unique_pops <- unique(mcols(pop_gr)$population_id)


# Sanitize function to make safe filenames
sanitize_filename <- function(name) {
  gsub("[^A-Za-z0-9_]", "_", name)  # Replace anything that's not a letter, number, or underscore
}

# Create directory for plots
dir.create("wild_roh_plots", showWarnings = FALSE)

dir.create("wild_roh_density_plots", showWarnings = FALSE) #skipped for now

# Loop through each sample and create sanitized output files
for (sample_id in unique_samples) {
  safe_id <- sanitize_filename(sample_id)
  output_file <- paste0("wild_roh_plots/", safe_id, "_roh_plot.pdf")
  plot_individual_roh(sample_id, output_file)
}

#loop through to create output files ~~~~### skipped for now
for (pop in unique_pops) {
  safe_name <- sanitize_filename(pop)
  output_file <- paste0("wild_roh_density_plots/pop_", safe_name, "_roh_density.pdf")
  plot_population_roh_smoothed(population_id = pop, output_file = output_file)
}
dev.off()

###ROH single chrom plot
#reading in data
roh_chr4 <- read.csv("hmm_roh_seg_categorized.csv", header = TRUE)
chr4_df <- data.frame(
  chr = paste0(roh_chr4$chromosome),  # Add 'chr' prefix if needed
  start = roh_chr4$start, #starting location of roh
  end = roh_chr4$end, #ending location of roh
  kb = roh_chr4$length_KB,   # ROH length in KB
  nsnp = roh_chr4$n_markers,   # Number of SNPs in ROH
  sample_id = roh_chr4$sample, # Individual ID
  pop = roh_chr4$pop #population information
)
#filtering dataframe for only chrome 10
chr4plot<- chr4_df[chr4_df$chr == "4", ]
#creating a custom genotype for karyoplotr for chrom 10 only
feline_chr_sizes <- data.frame(
  chr = "chr4", 
  start = 1,
  end = 207231548   # chr10 length only
) 
feline_genome_gr <- GRanges(seqnames = feline_chr_sizes$chr, ranges = IRanges(start = feline_chr_sizes$start, end = feline_chr_sizes$end)) #make into grange object
ocel_genome <- feline_genome_gr #making the custom plot type
seqlevels(ocel_genome) <- feline_chr_sizes$chr
seqlengths(ocel_genome) <- feline_chr_sizes$end
#making hom file into a granges file
chr4_gr <- GRanges(
  seqnames = paste0("chr", chr4plot$chr),
  ranges = IRanges(start = chr4plot$start, end = chr4plot$end),
  kb = chr4plot$kb,
  nsnp = chr4plot$nsnp,
  sample_id = chr4plot$sample_id)
#creating the function to plot
chr_4_plot <- function(sample_id, output_file = NULL) {
  # Filter ROH data for specific individual
  individual_roh <- chr4_gr[mcols(chr4_gr)$sample_id %in% sample_id]
  
  # Determine if output should go to a file
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 7)
  }
  #create the plot using the custom genome
  w_kp <- plotKaryotype(genome = ocel_genome, plot.type = 1, main = paste("ROH on Chromosome 4 for", sample_id))
  #kpAddChromosomeNames(w_kp, srt = 45, cex = 0.8) #not using this line for now
  #plot individuals
  kpRect(w_kp, 
         chr = as.character(seqnames(individual_roh)), 
         x0 = start(individual_roh), 
         x1 = end(individual_roh),
         y0 = 0, 
         y1 = 1, 
         col = "#FF000080",  # Semi-transparent red
         border = "#FF0000",
         r0 = 0.05, r1 = 0.95,
         data.panel = "ideogram") #should plot directly onto the ideogram
  # Close file if opened
  if (!is.null(output_file)) {
    dev.off()
  }
  # Return the filtered data
  return(individual_roh)
}
#using the function
unique_samples <- unique(chr4plot$sample_id)
sanitize_filename <- function(name) {
  gsub("[^A-Za-z0-9_]", "_", name)  # Replace anything that's not a letter, number, or underscore
}
dir.create("roh_plot_chr4", showWarnings = FALSE)
for (sample_id in unique_samples) {
  safe_id <- sanitize_filename(sample_id)
  output_file <- paste0("roh_plot_chr4/", safe_id, "_chr4_roh_plot.pdf")
  chr_4_plot(sample_id, output_file)
}




####Diversity stats -- Done!####

####expected, observed, and fis
##code from vcftools in terminal:
#vcftools --gzvcf wild_standard_final.vcf.gz --het  --out wild_het  
#read in output
ind_het <- read.table("wild_het.het", header = TRUE)
#converting raw values into the proportions
ind_het$O_het <- (ind_het$N_SITES - ind_het$O.HOM.) / ind_het$N_SITES
ind_het$O_hom <- ind_het$O.HOM. / ind_het$N_SITES
ind_het$E_hom <- ind_het$E.HOM. / ind_het$N_SITES
ind_het$E_het <- (ind_het$N_SITES - ind_het$E.HOM.) / ind_het$N_SITES
ind_het[, c("INDV", "O_hom", "O_het", "E_hom", "E_het", "F")]

write.csv(ind_het, "ind_het_fis.csv")

#get population averages
pop_id <- read.table("wild_pop_id.txt", header = FALSE)
colnames(pop_id) <- c("INDV", "pop")
as.data.frame(pop_id)
pop_het <- left_join(ind_het, pop_id, by = "INDV")

#mean expected het
pop_mean_he <- pop_het %>%
  group_by(pop) %>%
  summarise(Average = mean(E_het, na.rm = TRUE),
            Count = n(),
            StdDev = sd(E_het, na.rm = TRUE),
            StdError = (sd(E_het, na.rm = TRUE)/sqrt(n())))
write.csv(pop_mean_he, "pop_mean_he.csv")
#mean observed het
pop_mean_ho <- pop_het %>%
  group_by(pop) %>%
  summarise(Average = mean(O_het, na.rm = TRUE),
            Count = n(),
            StdDev = sd(O_het, na.rm = TRUE),
            StdError = (sd(O_het, na.rm = TRUE)/sqrt(n())))
write.csv(pop_mean_ho, "pop_mean_ho.csv")
#mean fis
pop_mean_f <- pop_het %>%
  group_by(pop) %>%
  summarise(Average = mean(F, na.rm = TRUE),
            Count = n(),
            StdDev = sd(F, na.rm = TRUE),
            StdError = (sd(F, na.rm = TRUE)/sqrt(n())))
write.csv(pop_mean_f, "pop_mean_F.csv")

###nucleotide diversity - population level
#command for producing windowed pi with vcftools, done for both populations:
#vcftools --gzvcf refuge_standard_postsubset.vcf.gz --window-pi 50000 --out refuge-nucleotide

#reading in output files
ranch_pi <- read.table("ranch-nucleotide.windowed.pi", header = TRUE)
refuge_pi <- read.table("refuge-nucleotide.windowed.pi", header = TRUE)
#adding population columns
ranch_pi$pop <- "Ranch"
refuge_pi$pop <- "Refuge"
#combine the two files
all_pi <- rbind(ranch_pi, refuge_pi)

#average pi per chromosome per population - weighted to reduce the affect of windows with low sites
chrom_pi <- all_pi %>%
  group_by(pop, CHROM) %>%
  summarise(mean_pi = weighted.mean(PI, N_VARIANTS, na.rm = TRUE),
            sd_pi = sd(PI, na.rm = TRUE),
            n_windows = n())
write_csv(chrom_pi, "chrom_pi_by_pop_weighted.csv")
#plot? trying to see if interesting at all
ggplot(chrom_pi, aes(x = factor(CHROM), y = mean_pi, fill = pop)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.00001),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Ranch" = "orchid", "Refuge" = "#01004c")) +
  labs(x = "Chromosome",
       y = "Mean π (Weighted)",
       fill = "Population") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 13.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))


#weighted genome wide pi -- weighted to reduce the affect of windows with low sites
genome_pi_weighted <- all_pi %>%
  group_by(pop) %>%
  summarise(mean_pi = weighted.mean(PI, N_VARIANTS, na.rm = TRUE),
            Count = n(),
            StdDev = sd(PI, na.rm = TRUE),
            StdError = (sd(PI, na.rm = TRUE)/sqrt(n())))
write.csv(genome_pi_weighted,"weighted_genome_pi.csv")


####coverage stats from bcftools - can delete later to clean####
###THIRD iteration, running diagnostics in bcftools

#getting distribution of coverage depths in bcftools:
  #cd /Users/tylerbostwick/bcftools #setting working directory for bcftools
  # ./bcftools stats /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/joint_call_autosomes.vcf.gz | grep "^SN"
        #gives summary stats on the vcf file
        #output:
        #SN	0	number of samples:	89
        #SN	0	number of records:	120040661 -- total number of sites
        #SN	0	number of no-ALTs:	0 -- rows with no alt allele, essentially uninformative rows
        #SN	0	number of SNPs:	101914847 -- total number of snps
        #SN	0	number of MNPs:	0 -- total of multi-nucleotide polys
        #SN	0	number of indels:	20246516 -- total number of indels (need to remove)
        #SN	0	number of others:	0 -- anything else
        #SN	0	number of multiallelic sites:	12342780 -- sites with more than two alleles
        #SN	0	number of multiallelic SNP sites:	1673221 -- number of multiallelic sites that are snps
#after removing all indels and biallelic sites -- total snps left is 103792041
#computing mean coverage depth across sites, code queries the vcf for the depth at each site, the performs a pipe that calculates mean
#./bcftools query -f '[%DP\n]' /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz | \ awk '{sum+=$1; n++} END {print "Mean per-individual depth:", sum/n}'
      #output: Mean per-individual depth: 12.0631
#computing mean and standard deviation bounds for use in determining the depth filter:
#./bcftools query -f '[%DP\n]' /Users/tylerbostwick/Documents/Masters_Work/Analyses/1_Data/1_Working_Files/Filter_testing/allInd_SNPs_autosomes.vcf.gz | \ awk '$1 > 0 {sum+=$1; sumsq+=$1*$1; n++} END {mean=sum/n; sd=sqrt(sumsq/n - mean^2); print "Mean:", mean; print "SD:", sd; print "Mean - 2SD:", mean-2*sd; print "Mean + 2SD:", mean+2*sd}'
      #output: Mean: 12.1187
              #SD: 5.08218
              #Mean - 2SD: 1.95431
              #Mean + 2SD: 22.2831






