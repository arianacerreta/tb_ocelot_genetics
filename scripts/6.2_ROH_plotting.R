####ROH -- Done!####
#set working directory to roh specific folder
setwd("./results/roh")

#clean your environment
rm(list = ls())

#libraries
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)
library(gridGraphics)

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
  ))

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
             color = "black", size = 4, shape = 18) +
  geom_errorbar(data = population_mean_roh, 
                aes(x = pop, y = Average, 
                    ymin = Average - 1*StdDev, ymax = Average + 1*StdDev), #updated to 1SD per comments-ALC9/9/2026
                color = "black", width = 0.2, size = 1) +
  scale_fill_manual(values = c("Ranch" = "#CC79A7", "Refuge" = "#009E73")) +
  # Labels and theme
  labs(x = "Population", y = expression(F[ROH] ("%"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none")

ggsave("../../figures/wild_roh_violin_plot.pdf", dpi = 1200)

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
sum(ind_roh_summary$Count)
mean(roh_seg_df$length_KB)

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
  scale_fill_manual(values = c("ranch" = "#CC79A7", "refuge" = "#009E73"),
                    labels = c("ranch" = "Ranch", "refuge" = "Refuge")) +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -10),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.margin = margin(b= 2, unit = "cm")
  )

ggdraw(p)+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[1],")"), x = 0.8, y = 0.13, size = 11, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[2],")"), x = 0.705, y = 0.13, size = 11, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[3],")"), x = 0.605, y = 0.13, size = 11, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[4],")"), x = 0.51, y = 0.13, size = 11, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[5],")"), x = 0.415, y = 0.13, size = 11, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[6],")"), x = 0.315, y = 0.13, size = 11, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[7],")"), x = 0.22, y = 0.13, size = 11, color = "black")+
  draw_label(paste0("(",(levels(unique(gen_plot_df$range_label)))[8],")"), x = 0.12, y = 0.13, size = 11, color = "black")+
  draw_label("(ROH length in Mb)", x = 0.465, y = 0.022, size = 13, color = "black")

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

#test
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
kpAxis(ind_kp, ymin = 0, ymax= 1, data.panel="ideogram")
kpAddBaseNumbers(w_kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1,
                 minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
kpPoints(ind_kp, 
         chr = as.character(seqnames(individual_roh)),
         x=chr4$window_mid,
         y=chr4$het_per_kb,
         data.panel= "ideogram"
)

