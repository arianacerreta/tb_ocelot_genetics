####ROH -- Done!####
#set working directory to roh specific folder
setwd("./results/roh")

#clean your environment
rm(list = ls())

library(data.table)

#libraries
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)
library(gridGraphics)
library(qpdf)
library(stringr)
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
###brief view of the output and froh by individual:
#read in output table, only the RG (roh segment) lines
roh_LD <- read.table("roh_LD05.txt", comment.char = "#", header = FALSE)
# Name the columns
colnames(roh_LD) <- c("type", "sample", "chromosome", "start", "end", "length_bp", "n_markers", "quality") #quality is average fwd-bwd phred score
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
ind_roh_LD$FROH <- (ind_roh_LD$length_bp / 2468705656) * 100 #updated number of bases to reflect ocelot scaffold primary assembly length (bp) Foley et al. 2026

#writing tables
write.csv(ind_roh_LD, "froh_by_individual_hmm.csv")
write.csv(roh_LD_filt, "roh_segments_filt_qual_length_hmm.csv", row.names = FALSE)

#adding population labels to data
pop_id<-read.delim("../../data/inputs/wild_subset.txt", header = FALSE, sep = "\t") #read in wild_subset.txt
pop_id<-pop_id[,-1]
colnames(pop_id)<- c("sample", "pop")
wild_roh_pop <- left_join(ind_roh_LD, pop_id, by = "sample")

#change pop labels to capitalized values for plotting
wild_roh_pop<-wild_roh_pop%>%
  mutate(pop = case_when(
    pop == "refuge" ~ "Refuge",
    pop == "ranch" ~ "Ranch",
    TRUE ~ pop
  ))%>%
  mutate(sample = gsub("-.*", "", sample))

#add in heterozygosity for each individual
#obs_het with same filtering as for Roh calculations 
pop_het_fis<-read.csv("../diversity/pop_het_fis_rohfilt.csv", header = TRUE)

wild_roh_pop<-wild_roh_pop%>%
  left_join(pop_het_fis, by = c("sample" = "INDV"))

FROHvHet<-ggplot(data = wild_roh_pop, aes(x = O_het, y = FROH, color = pop.x))+
  geom_point()+
  theme_minimal() +
  labs(x = "Observed heterozygosity", y = expression(F[ROH] ("%")), color = "Population")+
  scale_color_manual(values = c("Ranch" = "#CC79A7", "Refuge" = "#0072B2"))+
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = c(0.70,0.98),
        legend.justification = c("left", "top"),
        legend.box.background = element_rect(color = "#e5e5e5", fill = "#f9f9f9", linewidth = 0.5),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11),
        )

####ROH assessment
#average % genome in roh by population
population_mean_roh <- wild_roh_pop %>%
  group_by(pop.x) %>%
  summarise(Average = mean(FROH, na.rm = TRUE),
            Count = n(),
            StdDev = sd(FROH, na.rm = TRUE),
            StdError = (sd(FROH, na.rm = TRUE)/sqrt(n())))

write.csv(population_mean_roh, "population_mean_froh_table.csv")

#violin plot of wild FROH
ggplot() +
  # violin plot
  geom_violin(data = wild_roh_pop, aes(x = pop.x, y = FROH, fill = pop.x), 
              alpha = 0.7) +
  # Add population means with error bars
  geom_point(data = population_mean_roh, aes(x = pop.x, y = Average), 
             color = "black", size = 4, shape = 18) +
  geom_errorbar(data = population_mean_roh, 
                aes(x = pop.x, y = Average, 
                    ymin = Average - 1*StdDev, ymax = Average + 1*StdDev), #updated to 1SD per comments-ALC9/9/2026
                color = "black", width = 0.2, size = 1) +
  scale_fill_manual(values = c("Ranch" = "#CC79A7", "Refuge" = "#0072B2")) +
  # Add individual points
  geom_jitter(data = wild_roh_pop, aes(x = pop.x, y = FROH), #consider seeing if height = 0 changes anything
              width = 0.15, height = 0, alpha = 0.4, size = 3) +
  
  # Labels and theme
  labs(x = "Population", y = expression(F[ROH] ("%"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")

ggsave("../../figures/wild_roh_violin_plot.pdf", dpi = 1200) #change to .pdf for pdf; .png for png

####proportion of ROH lengths
#read in .hom files
roh_seg_df <- roh_LD_filt
#adding population information
roh_seg_df <- merge(roh_seg_df, pop_id, by = "sample", all = TRUE)
#calculate length in Mb -- bp to mb conversion
roh_seg_df$length_MB <- roh_seg_df$length_bp/1000000
roh_seg_df$length_KB <- roh_seg_df$length_bp/1000 #converting bp to kilobase pair
#look at resulting distribution
hist(roh_seg_df$length_MB, main="Distribution of ROH lengths", xlab="Length (MB)")
max(roh_seg_df$length_MB, na.rm = TRUE) #139.6876
min(roh_seg_df$length_MB, na.rm = TRUE) #0.526017
mean(roh_seg_df$length_MB, na.rm = TRUE) #11.19929
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
ind_roh_summary$FROH <- (ind_roh_summary$total_length / 2468705656) * 100 #updated number of bases to reflect ocelot scaffold primary assembly length (bp) Foley et al. 2026
sum(ind_roh_summary$Count) #2177
mean(roh_seg_df$length_KB) #11302.36

ind_bin_summary <- roh_seg_df %>% 
  summarize(.by = c(sample, Category),
    count = n()
  ) %>%
  arrange(sample, Category)

View(ind_bin_summary)

#plot generation time to roh segment
roh_seg_cat <- roh_seg_df
roh_seg_cat$gen <- 100/(2*(1.1)*roh_seg_cat$length_MB) #calculating generation time
#1.1 is recombination rate of the domestic cat (1.1 cM per Mb)

#binning generation time
roh_seg_cat <- roh_seg_cat %>% 
  mutate(
    gen_bins = case_when(
      gen < 2 ~ "<2",
      gen >= 2 & gen < 5 ~ "2-5",
      gen >= 5 & gen < 10 ~ "5-10",
      gen >= 10 & gen < 20 ~ "10-20",
      gen >= 20 & gen < 40 ~ "20-40",
      gen >= 40 & gen < 60 ~ "40-60",
      gen >= 60 & gen < 80 ~ "60-80",
      gen >= 12 ~ "80+",
      TRUE ~"Other"
    ),
    gen_bins = factor(gen_bins, levels = c("<2","2-5","5-10","10-20","20-40","40-60", "60-80","80+"))
  )

#get ranges of roh length in each generation bin; this is for data present
bin_ranges <- roh_seg_cat %>% 
  group_by(gen_bins) %>% 
  summarise(
    min_length = round(min(length_MB, na.rm =TRUE), 2),
    max_length = round(max(length_MB, na.rm = TRUE),2),
    range_span = max_length - min_length,
    .groups = 'drop'
  )

#back calculate hypothetical Mb from generation boundaries
bin_bound<-c(2,5,10,20,40,60,80)
MB<-round(100/(2*(1.1)*bin_bound),2)

bin_ranges$max_calc<-c(MB[1],MB)
bin_ranges$min_calc<-c(MB,MB[length(MB)])


#merge bin ranges with data frame
gen_plot_df <- roh_seg_cat %>% 
  left_join(bin_ranges, by = "gen_bins") %>% 
  mutate(range_label = paste0(min_calc, "-", max_calc))%>%
  mutate(range_label = case_when(
    range_label == "0.57-0.57" ~ "<0.57",
    range_label == "22.73-22.73" ~ ">22.73",
    TRUE ~ range_label
  ),
  range_label = factor(range_label, levels = c("<0.57","0.57-0.76","0.76-1.14","1.14-2.27","2.27-4.55","4.55-9.09", "9.09-22.73",">22.73"))
  )

#plotting generation time by population
p<-ggplot(gen_plot_df, aes(x = gen_bins, fill = pop)) +
  geom_bar(position = "dodge", stat = "count") +
  scale_y_continuous(breaks = breaks_width(100)) +
  labs(
    x = "Expected Generations",
    y = "ROH Count",
    fill = "Populations"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("ranch" = "#CC79A7", "refuge" = "#0072B2"),
                    labels = c("ranch" = "Ranch", "refuge" = "Refuge")) +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, vjust = -6),
    axis.title.y = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.box.background = element_rect(color = "#e5e5e5", fill = "#f9f9f9", size = 0.5),
    legend.position = c(0.98,0.98),
    legend.justification = c("right", "top"),
    plot.margin = margin(b= 1.5, t = 0.5,l = 0.7, unit = "cm")
  )

exp_gen<-ggdraw(p)+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[1],")"), x = 0.932, y = 0.23, size = 6, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[2],")"), x = 0.824, y = 0.23, size = 6, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[3],")"), x = 0.713, y = 0.23, size = 6, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[4],")"), x = 0.602, y = 0.23, size = 6, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[5],")"), x = 0.493, y = 0.23, size = 6, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[6],")"), x = 0.381, y = 0.23, size = 6, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[7],")"), x = 0.270, y = 0.23, size = 6, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[8],")"), x = 0.162, y = 0.23, size = 6, color = "black")+
  draw_label("(ROH length in Mb)", x = 0.55, y = 0.09, size = 8, color = "black")

 
#new Categories for Table 2
roh_seg_df$Category2 <- cut(roh_seg_df$length_MB,
                           breaks = c(0, 1, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70, 80, Inf),
                           labels = c("<1Mb", "1-2Mb", "2-4Mb", "4-6Mb", "6-8Mb", "8-10MB", "10-20Mb",
                                      "20Mb-30Mb", "30Mb-40Mb", "40Mb-50Mb", "50Mb-60Mb", "60Mb-70Mb",
                                      "70Mb-80Mb", ">80Mb"),
                         include.lowest = TRUE)
#summary by population  
population_category_counts <- roh_seg_df %>%
  group_by(pop, Category2)  %>%
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
    Avg_ROH_Count_Per_Individual = round((Count / n_individuals), digits = 2),
    Avg_ROH_Length_MB_Per_Individual = round((Total_Length_MB / n_individuals), digits = 2)
  ) %>%
  ungroup()
write.csv(normalized_roh, "pop_normalized_roh.csv")

##################################################################
####Visualizing ROH --karyotype plot --- individual plots working

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
  end = c(241490831,  # chr1 Chr lengths from NCBI ocelot reference (https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_058742865.1/)
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

##using the function -- plotting individual roh
unique_samples <- unique(plot_df$sample_id)

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
files<-list.files("./wild_roh_plots/")
files_w_path<-rep(paste0("./wild_roh_plots/",files))

pdf_combine(input = files_w_path, output = "./wild_roh_plots/all_indiv_roh_plots.pdf")

#clean enviornment of all but the plots needed to combine at end
rm(list=setdiff(ls(),c("FROHvHet", "exp_gen")))

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
#filtering dataframe for only chrome 4
chr4plot<- chr4_df[chr4_df$chr == "4", ]

#creating a custom genotype for karyoplotr for chrom 4 only
feline_chr_sizes <- data.frame(
  chr = "chr4", 
  start = 1,
  end = 207231548  # chr4 length only 

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
files<-list.files("./roh_plot_chr4/")
files_w_path<-rep(paste0("./roh_plot_chr4/",files))

pdf_combine(input = files_w_path, output = "./roh_plot_chr4/all_indiv_chr4_roh_plots.pdf")

#plotting for Figure 3
#creating a custom genotype for karyoploter for chrom 4 only
#I have to duplicate chr4 for the total number of individuals I want to
feline_chr_sizes <- data.frame(
  chr = c("LO01F", "LO03M", "OM283", "OF304", "E14F", "E33M"),
  start = rep(1,6),
  end = rep(207231548,6)  # chr4 length only
) 
feline_genome_gr <- GRanges(seqnames = feline_chr_sizes$chr, ranges = IRanges(start = feline_chr_sizes$start, end = feline_chr_sizes$end)) #make into grange object
ocel_genome <- feline_genome_gr #making the custom plot type
seqlevels(ocel_genome) <- feline_chr_sizes$chr
seqlengths(ocel_genome) <- feline_chr_sizes$end

#making hom file into a granges file
chr4plot<-chr4plot %>%
  mutate(sample_id = gsub("-.*", "", sample_id))
chr4_gr <- GRanges(
  seqnames = paste0(chr4plot$sample_id),
  ranges = IRanges(start = chr4plot$start, end = chr4plot$end),
  kb = chr4plot$kb,
  nsnp = chr4plot$nsnp,
  sample_id = chr4plot$sample_id)

#######
#calculate windowed heterozygosity (code pulled from Matt Smith)
# Set parameters
window_size <- 1e6             # 1 Mb
step_size <- window_size       # non-overlapping
low_het_threshold <- 0.005     # cutoff for low het

# Define chromosome order (adjust to match your labeling scheme if needed)
chr_order <- c(1:17)

# Read in PLINK raw file (genotype: 0=homo1, 1=het, 2=homo2)
geno <- fread("../../data/processed/addtl_filter/roh_LDpruned_05_data_allele.raw", header = TRUE) # Use fread to handle large tables
geno_data <- geno %>% select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

# Extract SNP metadata (assumes format like chr:pos for SNP names)
snp_info <- data.frame(snp = colnames(geno_data)[-1], stringsAsFactors = FALSE) %>%
  mutate(chr = str_split_i(snp, "_", 1), pos = as.integer(str_split_i(snp, "_",2)))

# Ensure chromosome ordering is preserved
snp_info$chr <- factor(snp_info$chr, levels = chr_order)

# Genotype matrix to long format
geno_long <- geno_data %>%
  dplyr::rename(individual = IID) %>%
  tidyr::pivot_longer(
    cols = -individual,
    names_to = "snp",
    values_to = "geno"
  ) %>%
  left_join(snp_info, by = "snp")

# Remove large dataframes
rm(geno); rm(geno_data)

# Filter out SNPs with no position info
geno_long <- geno_long %>% filter(!is.na(pos))

#filter to chr4 to figure out what is going on
geno_long_chr4<- geno_long %>% filter(chr =="4")
max(geno_long_chr4$pos) #207,193,806

# Calculate windowed heterozygosity by individual and chromosome
windowed_het <- data.frame()

for (indiv in unique(geno_long$individual)) {
  indiv_data <- geno_long %>% filter(individual == indiv)
  
  for (chr in levels(snp_info$chr)) {
    chr_data <- indiv_data %>% filter(chr == !!chr)
    if (nrow(chr_data) == 0) next
    
    max_pos <- max(chr_data$pos, na.rm = TRUE)
    starts <- seq(0, max_pos, by = step_size)
    
    for (start in starts) {
      end <- start + window_size
      window <- chr_data %>% filter(pos >= start & pos < end)
      window_geno <- window$geno
      n_het <- sum(window_geno == 1, na.rm = TRUE)
      n_snps <- sum(!is.na(window_geno))
      het_per_kb <- n_het / (window_size / 1000)
      het_per_snp <- ifelse(n_snps > 0, n_het / n_snps, NA)
      
      windowed_het <- rbind(windowed_het, data.frame(
        individual = indiv,
        chr = chr,
        window_start = start,
        window_end = end,
        window_mid = start + window_size / 2,
        n_het = n_het,
        n_snps = n_snps,
        het_per_kb = het_per_kb,
        het_per_snp = het_per_snp
      ))
    }
  }
}

# Set chromosome as factor with defined order
windowed_het$chr <- factor(windowed_het$chr, levels = chr_order)

#clean individual names
windowed_het<-windowed_het %>%
  mutate(individual = gsub("-.*", "", individual))

# Compute average across individuals
windowed_het_avg <- windowed_het %>%
  group_by(chr, window_mid) %>%
  summarise(
    mean_het_per_kb = mean(het_per_kb, na.rm = TRUE),
    mean_het_per_snp = mean(het_per_snp, na.rm = TRUE),
    .groups = "drop"
  )

# Add cumulative genome position for plotting
chr_info <- windowed_het %>%
  group_by(chr) %>%
  summarise(chr_len = max(window_end), .groups = "drop") %>%
  mutate(chr_start = cumsum(lag(chr_len, default = 0)))

windowed_het <- windowed_het %>%
  left_join(chr_info, by = "chr") %>%
  mutate(cum_pos = window_mid + chr_start)

windowed_het_avg <- windowed_het_avg %>%
  left_join(chr_info, by = "chr") %>%
  mutate(cum_pos = window_mid + chr_start)

chr_midpoints <- chr_info %>% mutate(mid = chr_start + chr_len / 2)

#filter to only Chr4
windowed_het_chr4<-windowed_het %>%
  filter(chr == "4")
LO01F<-windowed_het_chr4 %>%
  filter(individual == "LO01F")

saveRDS(windowed_het, file = "windowed_het.RDS")

#read in RDS if needed
windowed_het<-readRDS(file = "windowed_het.RDS")
#make a GRange object
het_points <- GRanges(
  seqnames = paste0(windowed_het_chr4$individual), #since this is all chr4; hacking this line for individual labels
  ranges = IRanges(start = windowed_het_chr4$window_mid, width = rep(1,length(windowed_het_chr4$window_mid))),
  het= windowed_het_chr4$het_per_snp,#heterozygous sites in a window/called SNP positions in a window; ala Saremi et al 2019 Supplemental Methods
  sample_id = windowed_het_chr4$individual)

#LO01F-1,LO03M-1,OM283-2 (big ROH),LAO03M-1 (big ROH), OF304-1A (little ROH), E33M-1 (little ROH), E14F-1 (big ROH)
sample_id<-c("LO01F","LO03M", "OM283", "OF304", "E14F", "E33M") #ranch
individual_roh <- chr4_gr[mcols(chr4_gr)$sample_id %in% sample_id]
#"ranch" = "#CC79A7", "refuge" = "#0072B2"
color_ref<-matrix(data = NA, ncol = 3, nrow = length(individual_roh@seqnames@values))
colnames(color_ref)<-c("sample_id", "seg", "color")
color_ref<-as.data.frame(color_ref)
color_ref$sample_id<-paste0(individual_roh@seqnames@values)
color_ref$seg<-individual_roh@seqnames@lengths
color_ref$color<-c("#CC79A7","#CC79A7","#0072B2", "#0072B2", "#0072B2","#0072B2")
colors<-rep(color_ref$color, color_ref$seg)

params<-getDefaultPlotParams(6)
params$topmargin<-1
params$bottommargin<-3
params$ideogramheight<-5
params$dataideogrammax<-1
params$data1outmargin<-0
params$data2outmargin<-2
params$leftmargin<-0.15

#start plot
w_kp <- plotKaryotype(genome = ocel_genome, chromosomes = "all",
                      plot.params = params,
                      plot.type = 6,
                      cex = 0.5
                      )
kpDataBackground(w_kp,data.panel = "ideogram", color = "white", r0=0.0, r1=0.98, clipping = FALSE)

kpRect(w_kp, 
       chr = as.character(seqnames(individual_roh)), 
       x0 = start(individual_roh), 
       x1 = end(individual_roh),
       y0 = 0, 
       y1 = 1, 
       col = colors, 
       border = "black",
       r0 = 0.00, r1 = 1,
       data.panel = "ideogram") #should plot directly onto the ideogram
kpAddBaseNumbers(w_kp, tick.dist = 100000000, tick.len = 0.5, tick.col="black", cex=0.5,
                 minor.tick.dist = 10000000, minor.tick.len = 0, minor.tick.col = "black")

kpAxis(w_kp, ymin = 0, ymax= 1, cex = 0.6,data.panel="ideogram")

kpPoints(w_kp, 
         data=het_points,
         y=het_points$het,
         data.panel= "ideogram",
         cex =0.4
)
#end plot
kary<-recordPlot()

kary_labelled<-ggdraw(kary)+
  draw_label("Position on Chr 4 (Mb)", x = 0.55, y = 0.05, size = 11, color = "black")+
  draw_line(x=c(0.075,0.075), y=c(0,1), color = "white", size = 20)+
  draw_line(x=c(0,1), y=c(0.24,0.24), color = "white", size = 3)+
  draw_line(x=c(0,1), y=c(0.39,0.39), color = "white", size = 3)+
  draw_line(x=c(0,1), y=c(0.54,0.54), color = "white", size = 3)+
  draw_line(x=c(0,1), y=c(0.69,0.69), color = "white", size = 3)+
  draw_line(x=c(0,1), y=c(0.84,0.84), color = "white", size = 3)+
  draw_label("1", x=0.135, y=0.976, color = "black", size = 6)+
  draw_label("1", x=0.135, y=0.826, color = "black", size = 6)+
  draw_label("1", x=0.135, y=0.676, color = "black", size = 6)+
  draw_label("1", x=0.135, y=0.526, color = "black", size = 6)+
  draw_label("1", x=0.135, y=0.374, color = "black", size = 6)+
  draw_label("1", x=0.135, y=0.224, color = "black", size = 6)+
  draw_label("0.5", x=0.127, y=0.922, color = "black", size = 6)+
  draw_label("0.5", x=0.127, y=0.772, color = "black", size = 6)+
  draw_label("0.5", x=0.127, y=0.622, color = "black", size = 6)+
  draw_label("0.5", x=0.127, y=0.472, color = "black", size = 6)+
  draw_label("0.5", x=0.127, y=0.32, color = "black", size = 6)+
  draw_label("0.5", x=0.127, y=0.17, color = "black", size = 6)+
  draw_label("0", x=0.135, y=0.868, color = "black", size = 6)+
  draw_label("0", x=0.135, y=0.718, color = "black", size = 6)+
  draw_label("0", x=0.135, y=0.566, color = "black", size = 6)+
  draw_label("0", x=0.135, y=0.415, color = "black", size = 6)+
  draw_label("0", x=0.135, y=0.265, color = "black", size = 6)+
  draw_label("0", x=0.135, y=0.114, color = "black", size = 6)+
  draw_label("Heterozygosity", angle = 90, x = 0.09, y=0.55, size = 11, color = "black")+
  draw_label("LO01F", x=0.535, y=0.965, color = "black", size = 6)+
  draw_label("LO03M", x=0.535, y=0.81, color = "black", size = 6)+
  draw_label("OM283", x=0.535, y=0.66, color = "black", size = 6)+
  draw_label("OF304", x=0.535, y=0.51, color = "black", size = 6)+
  draw_label("E14F", x=0.536, y=0.357, color = "black", size = 6)+
  draw_label("E33M", x=0.536, y=0.207, color = "black", size = 6)


g1<-plot_grid(kary_labelled, FROHvHet,
              ncol = 2,
              labels = c('A','B'),
              rel_widths = c(1.25,1))

figure2<-plot_grid(g1,exp_gen,
                   nrow = 2,
                   labels = c("","C"), 
                   rel_heights = c(1.25,1))

ggsave("../../figures/figure3.pdf", dpi = 1200)
ggsave("../../figures/figure3.png", dpi = 1200)
