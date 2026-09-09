#Kinship analysis
#library
library(dplyr)
library(ggplot2)

#clean your environment
rm(list = ls())

#working directory
setwd("./data/processed/addtl_filter/")

#Define paths
##update with your path to PLINK and PLINK2
PLINKpath<-"F:/2_TB_Working_Files/Plink_files/WindowsPLINK"

##king-robust kinship estimator for WILD individuals 
system(paste0(PLINKpath,"/plink2 --bfile wild_kin_roh_filter --make-king-table --allow-extra-chr --chr-set 17 --out wild_KING_manu"))

king.wild.matrix <- read.table("wild_KING_manu.kin0", header = FALSE) #make king table into object
colnames(king.wild.matrix) <- c("FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") #add column header

wild_origins<-read.delim("../../inputs/wild_subset.txt", header = FALSE, sep = "\t") #read in wild_subset.txt
wild_origins<-wild_origins[,-1]
colnames(wild_origins)<- c("ID", "Pop")
wild_origins <- wild_origins %>%
  mutate(ID = gsub("-.*", "", ID))

#remove "-x" from the individual id's
king.wild.matrix <- king.wild.matrix %>%
  mutate(IID1 = gsub("-.*", "", IID1),
         IID2 = gsub("-.*", "", IID2))
king.wild.matrix <- left_join(king.wild.matrix, wild_origins, by = c("IID1"="ID"))%>%
  rename(Pop1=Pop)

king.wild.matrix<-left_join(king.wild.matrix, wild_origins, by = c("IID2"="ID"))%>%
  rename(Pop2=Pop)

king.wild.matrix<-king.wild.matrix%>%
  mutate(Pop1 = case_when(
    IID1 == "OM331" ~ "rdisp",
    IID1 == "LO03M" ~ "rdisp",
    IID1 == "E29M" ~ "rdisp",
    IID1 == "E32M" ~ "rdisp",
    TRUE ~ Pop1
  ))%>%
  mutate(Pop2 = case_when(
    IID2 == "OM331" ~ "rdisp",
    IID2 == "LO03M" ~ "rdisp",
    IID2 == "E29M" ~ "rdisp",
    IID1 == "E32M" ~ "rdisp",
    TRUE ~ Pop2
  ))

king.wild.matrix1<-king.wild.matrix%>%
  arrange(Pop1)

ordered_ID<-c("E10F",unique(king.wild.matrix1$IID1))

ordered_IID1<-ordered_ID[2:length(ordered_ID)]

ordered_IID2<- ordered_ID[1:length(ordered_ID)-1]

#make a new matrix so to organize the kinship values and IID1 and IID2 according to ranch, disp, refuge for manuscript
#this is becuase some IID1, IID2 pairs need to be swapped for proper plotting
new_ordered_matrix_king<-as.data.frame(matrix(data = NA, nrow=length(king.wild.matrix$IID1), ncol = 3))
colnames(new_ordered_matrix_king)<-c("IID1", "IID2", "KINSHIP")

IID1<-rep(ordered_IID1, times=1:length(ordered_IID1))
IID2<-vector()
for (i in 1:length(ordered_IID2)) {
  IDs<-ordered_IID2[1:i]
  IID2<-append(IID2,IDs)
}

new_ordered_matrix_king$IID1<-IID1
new_ordered_matrix_king$IID2<-IID2

for (i in 1:length(new_ordered_matrix_king$IID1)){
  IDs<-c(new_ordered_matrix_king[i,1], new_ordered_matrix_king[i,2])
  Kin<-king.wild.matrix%>%
  filter((IID1 %in% IDs[1]) & (IID2 %in% IDs[2]) |
           (IID2 %in% IDs[1]) & (IID1 %in% IDs[2]))
  new_ordered_matrix_king[i,3]<-Kin$KINSHIP
}

new_ordered_matrix_king<- new_ordered_matrix_king %>%
  mutate(IID1 = factor(IID1, levels = ordered_IID1))%>%
  mutate(IID2 = factor(IID2, levels = ordered_IID2))

#make csv from the .kin0 file
dir.create("../../../results/kinship")
write.csv(king.wild.matrix, "../../../results/kinship/wild_pairwise_kinship_manu.csv")

#heat map of kinship for wild individuals -- binned
ggplot(data = new_ordered_matrix_king, aes(x=IID1, y=IID2, fill = KINSHIP)) +
  geom_tile(color="white") +
  scale_fill_stepsn(name = "Kinship", breaks = c(0, 0.04, 0.1, 0.2),
                    limit = c(0, 0.3),
                    labels = c("Unrelated", "3rd Degree", "2nd Degree", "1st Degree"),
                    aesthetics = "fill",
                    space = "Lab",
                    colors = c("grey100", "gold1", "orangered", "firebrick"),
                    )+
  annotate("rect", xmin=c(29.5,29.5), xmax=c(30.5,30.5), ymin=c(13.5, 25.5), ymax=c(14.5,26.5),
          color = "black", fill = "transparent", size =1.5, linetype="solid" )+
  annotate("rect", xmin=c(0.5), xmax=c(22.5), ymin=c(0.5), ymax=c(23.5),
           color = "black", fill = "transparent", size =0.75, linetype="solid" )+
  annotate("text", x = 8, y = 18, label = "Ranch" , size = 5)+
  annotate("rect", xmin=c(26.5), xmax=c(43.5), ymin=c(27.5), ymax=c(43.5),
           color = "black", fill = "transparent", size =0.75, linetype="solid" )+
  annotate("text", x = 33, y = 41, label = "Refuge" , size = 5)+
  annotate("rect", xmin=c(22.5), xmax=c(26.5), ymin=c(23.5), ymax=c(27.5),
           color = "black", fill = "transparent", size =0.75, linetype="solid" )+
  annotate("text", x = 19.5, y = 28, label = "Dispersers" , size = 5)+
  labs(x="Individuals", y="Individuals")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(hjust = 1.25),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
