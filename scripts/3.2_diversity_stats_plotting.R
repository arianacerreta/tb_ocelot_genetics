#Diversity Stats and Coverage Stats
####expected, observed, and fis
library(ggplot2)
#use output from 3.1_bcftols_het.sh
#read in output
ind_het <- read.table("./results/diversity/wild_het.het", header = TRUE)
#converting raw values into the proportions

ind_het$O_het <- (ind_het$N_SITES - ind_het$O.HOM.) / ind_het$N_SITES
ind_het$O_hom <- ind_het$O.HOM. / ind_het$N_SITES
ind_het$E_hom <- ind_het$E.HOM. / ind_het$N_SITES
ind_het$E_het <- (ind_het$N_SITES - ind_het$E.HOM.) / ind_het$N_SITES
ind_het[, c("INDV", "O_hom", "O_het", "E_hom", "E_het", "F")]

write.csv(ind_het, "./results/diversity/ind_het_fis.csv")

#get population averages
wild_origins<-read.delim("./data/inputs/wild_subset.txt", header = FALSE, sep = "\t") #read in wild_subset.txt
wild_origins<-wild_origins[,-1]
colnames(wild_origins) <- c("INDV", "pop")
pop_het <- left_join(ind_het, wild_origins, by = "INDV")
pop_het <- pop_het %>%
  mutate(INDV = gsub("-.*", "", INDV))
write.csv(pop_het, "./results/diversity/pop_het_fis.csv")
#mean expected het
pop_mean_he <- pop_het %>%
  group_by(pop) %>%
  summarise(Average = mean(E_het, na.rm = TRUE),
            Count = n(),
            StdDev = sd(E_het, na.rm = TRUE),
            StdError = (sd(E_het, na.rm = TRUE)/sqrt(n())))
write.csv(pop_mean_he, "./results/diversity/pop_mean_he.csv")
#mean observed het
pop_mean_ho <- pop_het %>%
  group_by(pop) %>%
  summarise(Average = mean(O_het, na.rm = TRUE),
            Count = n(),
            StdDev = sd(O_het, na.rm = TRUE),
            StdError = (sd(O_het, na.rm = TRUE)/sqrt(n())))
write.csv(pop_mean_ho, "pop_mean_ho.csv")
#mean fis (based on Allendorf et al. 2022 this is actually just F=1-(Ho/He))
pop_mean_f <- pop_het %>%
  group_by(pop) %>%
  summarise(Average = mean(F, na.rm = TRUE),
            Count = n(),
            StdDev = sd(F, na.rm = TRUE),
            StdError = (sd(F, na.rm = TRUE)/sqrt(n())))
write.csv(pop_mean_f, "pop_mean_F.csv")
#Fis isn't calculated by averaging F which is 1-(Ho/He)
#Fis = 1- (Ho/Hs); Ho is obs het averaged over all subpopulations (or in this case if we are looking at ranch v refuge then one or the other)
#Hs is expected het averaged over all populations
#I also don't think we report an SE typically
overall_FIS<- 1-(mean(pop_het$O_het)/mean(pop_het$E_het))
ranch_FIS<-1-(pop_mean_ho[1,2]/pop_mean_he[1,2])
refuge_FIS<- 1-(pop_mean_ho[2,2]/pop_mean_he[2,2])


###nucleotide diversity - population level
#reading in output files
ranch_pi <- read.table("./results/diversity/ranch-nucleotide.windowed.pi", header = TRUE)
refuge_pi <- read.table("./results/diversity/refuge-nucleotide.windowed.pi", header = TRUE)
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
write.csv(chrom_pi, "./results/diversity/chrom_pi_by_pop_weighted.csv")

#plot
ggplot(chrom_pi, aes(x = factor(CHROM), y = mean_pi, fill = pop)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.00001),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Ranch" = "#CC79A7", "Refuge" = "#0072B2")) +
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

ggsave("./figures/weighted_pi.png", dpi = 1200)

#looking at it with SD; don't think these are that different tbh
ggplot() +
  # Add population means with error bars
  geom_point(data = chrom_pi, aes(x = factor(CHROM), y = mean_pi, color = pop), 
             size = 4, shape = 18) +
  geom_errorbar(data = chrom_pi, 
                aes(x = factor(CHROM), y = mean_pi, color = pop,
                    ymin = mean_pi - 1*sd_pi, ymax = mean_pi + 1*sd_pi), #updated to 1SD per comments-ALC9/9/2026
                width = 0.2, size = 1) +
  scale_color_manual(values = c("Ranch" = "#CC79A7", "Refuge" = "#009E73")) 

#weighted genome wide pi -- weighted to reduce the affect of windows with low sites
genome_pi_weighted <- all_pi %>%
  group_by(pop) %>%
  summarise(mean_pi = weighted.mean(PI, N_VARIANTS, na.rm = TRUE),
            Count = n(),
            StdDev = sd(PI, na.rm = TRUE),
            StdError = (sd(PI, na.rm = TRUE)/sqrt(n())))
write.csv(genome_pi_weighted,"./results/diversity/weighted_genome_pi.csv")
