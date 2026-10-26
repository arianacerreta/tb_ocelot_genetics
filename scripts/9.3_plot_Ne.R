library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

###plotting Ne
setwd("./results/current_Ne/")

###compare maf
#Since GONE2 uses only polymorphic loci, pre-filtered with --mac 1 in plink
no_maf<-read.delim("./gone2_r1_1_x_E_M0_wild_kin_roh_filter_mac1_GONE2_mix_density.tsv")
# Number of markers in the input file: 477244
# Number of SNPs included in the analysis (only polymorphic and with less than 20% missing data):20000
maf01<-read.delim("./gone2_r1_1_x_E__M01_wild_kin_roh_filter_mac1_GONE2_mix_density.tsv")
# Number of markers in the input file: 477244
# Number of SNPs included in the analysis (only polymorphic and with less than 20% missing data): 179006
###don't understand why GONE2 didn't subsample to 20000 by default...I tried forcing with -s 20000 and it did the same thing
maf05<-read.delim("./gone2_r1_1_x_E__M05_wild_kin_roh_filter_mac1_GONE2_mix_density.tsv")
# Number of markers in the input file: 477244
# Number of SNPs included in the analysis (only polymorphic and with less than 20% missing data): 165855
###don't understand why GONE2 didn't subsample to 20000 by default...did not bother attempting to force since that didn't work for maf01
no_maf_NT<-read.delim("./gone2_r1_1_x_E_M0_wild_kin_roh_filter_mac1_GONE2_mix_Ne")
maf01_NT<-read.delim("./gone2_r1_1_x_E__M01_wild_kin_roh_filter_mac1_GONE2_mix_Ne")
maf05_NT<-read.delim("./gone2_r1_1_x_E__M05_wild_kin_roh_filter_mac1_GONE2_mix_Ne")

no_maf_NT$maf<-rep(0, length(no_maf_NT$Generation))
maf01_NT$maf<-rep(0.01, length(maf01_NT$Generation))
maf05_NT$maf<-rep(0.05, length(maf05_NT$Generation))


compare_maf<-rbind(no_maf_NT,maf01_NT, maf05_NT)
gens<-compare_maf$Generation
years<-7*gens
compare_maf$years<-years
compare_maf$maf<-as.character(compare_maf$maf)

no_maf$maf<-rep(0, length(no_maf$Generation))
maf01$maf<-rep(0.01,length(maf01$Generation))
maf05$maf<-rep(0.05, length(maf05$Generation))

combo_density<-rbind(no_maf,maf01,maf05)
combo_density$maf<-as.character(combo_density$maf)

summary_density<-combo_density%>%
  group_by(maf, Generation)%>%
  summarize(
    Count = n(),
    mean = mean(Ne),
    sd= sd(Ne),
    median = median(Ne),
    max = max(Ne),
    min = min(Ne)
  )

summary_all<- left_join(summary_density, compare_maf, by = c("maf","Generation"))

ggplot(summary_all, aes(x = years, y = N_T, color = maf))+
  geom_ribbon(aes(x = years, ymax = N_T + sd, ymin = N_T -sd, fill = maf), alpha = 0.3)+
  geom_step(linewidth = 1)+
  labs(x= "Years Ago", y = expression(N[e]), fill = "MAF", color = "MAF")+
  scale_color_manual(values = c("magenta","blue","yellow"),
                     limits = c("0", "0.01", "0.05"))+
  scale_fill_manual(values = c("magenta","blue","yellow"),
                    limits = c("0", "0.01", "0.05"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,400))+
  scale_x_continuous(expand = c(0,0))+
  theme_bw()+
  theme(
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11))

ggsave("../../figures/ne_for_maf.png", dpi = 1200) #change to .pdf for pdf; .png for png

#########################################


only_no_maf<- summary_all%>%
  filter(maf == "0")

Ne<-ggplot(only_no_maf, aes(x = years, y = N_T))+
  annotate("rect",xmin = (2024-1905), xmax = (2024-1853), ymin = 0, ymax =320, fill = "gray70",  alpha = 0.4 )+
  annotate("rect",xmin = (2024-1981), xmax = (2024-1905), ymin = 0, ymax = 320, fill = "#faa19c",  alpha = 0.4 )+
  annotate("rect", xmin = (2024-1492), xmax = 1050, ymin = 0, ymax = 320, fill = "gray70",  alpha = 0.4)+
  geom_ribbon(aes(x = years, ymax = N_T + sd, ymin = N_T -sd),  fill = "gray50")+
  geom_step(linewidth = 0.7, color = "black")+
  labs(x= "Years Ago", y = expression(N[e]))+
  scale_y_continuous(expand = c(0,0), limits = c(0,320))+
  scale_x_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11))

Ne_final<-ggdraw(Ne)+
  draw_label("Arkansas specimen collected", x = 0.25, y = 0.63, size = 7, color = "gray30", angle = 90)+
  draw_label("Suspected extirpation", x = 0.205, y = 0.63, size = 7, color = "#ca4239", angle = 90)+
  #draw_label("1492 CE", x = 0.53, y = 0.58, size = 10, color = "gray30", angle = 90)+
  draw_label("Pre-Columbian era", x = 0.58, y = 0.63, size = 7, color = "gray30", angle= 90)


