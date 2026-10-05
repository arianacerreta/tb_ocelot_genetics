##admixture plotting
#using LD pruned dataset
#.q files from ADMIXTURE performed on a mac computer/server

#packages
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(dplyr)
library(cowplot)

#rename and zip files for CLUMPAK (https://clumpak.evolseq.net/)
#first rename all the files in each run# folder
runs<-10 #put how many runs you did
setwd("./results/admixture")
for (i in 1:runs){
  setwd(paste0("./run",i))
  Q_files<-list.files(pattern = "\\.Q$")
  split_names<-str_split(Q_files, pattern = "[.]")
  new_Q_names<-vector()
  for (k in 1:length(split_names)){
    front<-split_names[[k]][1]
    k<-paste0("_k",split_names[[k]][2])
    run<-paste0("_run",i)
    new_name<-paste0(front,k,run,".Q")
    new_Q_names<-append(new_Q_names,new_name)
  }
  file.rename(from = Q_files, to = new_Q_names)
  setwd("..")
}

dir.create("./forCLUMPAK")
admix_dir<-getwd()
for (i in 1:runs){
  setwd(paste0("./run",i))
  Q_files<-list.files(pattern = "\\.Q$")
  current_dir<-rep(getwd(),length(Q_files))
  new_dir<-rep(paste0(admix_dir,"/forCLUMPAK"),length(Q_files))
  from<- paste0(current_dir,"/",Q_files)
  to<- paste0(new_dir,"/",Q_files)
  file.copy(from = from, to = to)
  setwd("..")
}


setwd("./forCLUMPAK")
zip("forCLUMPAK.zip", files = list.files())
setwd("../../..")

wild_origins<-read.delim("./data/inputs/wild_subset.txt", header = FALSE, sep = "\t") #read in wild_subset.txt
populations_file<-wild_origins[,-c(1:2)]

write.table(populations_file, file= "./results/admixture/forCLUMPAK/populations_file.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# use forCLUMPAK.zip and populations_file.txt at https://clumpak.evolseq.net/ to evaluate multiple modes for each K

#clean your environment
rm(list=setdiff(ls(),c("wild_pca", "kinship_plot"))) #all but pca, kinship plot

###wild admixture plot
#look at cross-validation errors for multiple runs
#pull .txt files
runs<-10 #put how many runs you did
setwd("./results/admixture")
cv_for_ks<-matrix(ncol = 3)
for (i in 1:runs){
  setwd(paste0("./run",i))
  file<-list.files(pattern = "^CV_per_K_run[0-9]+\\.txt$")
  temp<-read.delim(file, header = FALSE, sep = "")
  cv_for_ks<-rbind(cv_for_ks, temp)
  setwd("..")
}
cv_for_ks<-cv_for_ks[-1,]
colnames(cv_for_ks)<-c("K","CV","run")

cv_summary<-cv_for_ks %>%
  group_by(K) %>%
  summarise(Average = mean(CV, na.rm = TRUE),
            Count = n(),
            StdDev = sd(CV, na.rm = TRUE),
            StdError = (sd(CV, na.rm = TRUE)/sqrt(n())))


CV_plot<-ggplot()+
  geom_jitter(data = cv_for_ks, aes(x=K,y=CV), alpha = 0.3, height = 0, width = 0.1, size = 3)+
  geom_point(data = cv_summary, aes(x=K,y= Average),
             color = "#0072B2", size = 3)+
  geom_line(data = cv_summary, aes(x=K,y= Average),
            color = "#0072B2", size = 1)+
  theme_minimal() +
  labs(y="Cross Validation Error")+
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
CV_plot
ggsave("../../figures/CV_ADMIXTURE_PLOT.pdf", dpi = 1200)

#reading in and preparing data
w_k2_table <- read.table("./run1/wild_LDpruned_05_k2_run1.Q") ## finished through here

#add in population information
wild_origins<-read.delim("../../data/inputs/wild_subset.txt", header = FALSE, sep = "\t") #read in wild_subset.txt
wild_origins<-wild_origins[,-1]
colnames(wild_origins)<- c("ID", "Pop")
wild_origins<-wild_origins %>%
  mutate(ID = gsub("-.*", "", ID))
w_k2 <- cbind(wild_origins, w_k2_table)
w_k2_long <- pivot_longer(
  w_k2, cols = c(V1, V2),
  names_to = "ancestry",
  values_to = "proportion"
) #pivots the table to the proportions are able to be plotted as stacked bars

w_k2_long$Pop <- str_to_title(w_k2_long$Pop)

#plot code
w_k2plot <-
  ggplot(w_k2_long, aes(x = factor(ID), y = proportion, fill = ancestry)) +
  geom_col(color = "gray", linewidth = 0.1) +
  facet_grid(~fct_inorder(Pop), switch = "both", scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, size = 12),
        panel.grid = element_blank(),
     #   axis.ticks.y.left = element_line(linewidth = 0.75, color = "gray"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        strip.text.x = element_text(face = "bold", size = 12, vjust = 2),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8, hjust = -1)
  ) +
  scale_fill_manual(values = c("V1" = "#0072B2", "V2" = "#CC79A7")) +
  guides(fill = "none")


w_k2plot

ggsave("../../figures/wild__k2_admixture_plot.pdf", dpi = 1200)

#clean environment
rm(list=setdiff(ls(),c("wild_pca", "kinship_plot", "w_k2plot"))) #all but plots

g1<-plot_grid(wild_pca, kinship_plot,
          ncol = 2,
          labels = c('A','B'),
          rel_widths = c(1,1))

figure2<-plot_grid(g1,w_k2plot,
          nrow = 2,
          labels = c("","C"), 
          rel_heights = c(2,1))

ggsave("../../figures/figure2.pdf", dpi = 1200)
ggsave("../../figures/figure2.png", dpi = 1200)
