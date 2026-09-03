#PCA Analysis
##PCA on Wild individuals

##Set paths
PLINKpath<-"F:/2_TB_Working_Files/Plink_files/WindowsPLINK"
BFILEpath<-"./data/processed/wild_standard_final"
OUTPUTpath<-"./results/pca"

###reading in data
system(paste0(PLINKpath,"/plink --bfile ", BFILEpath," --pca --chr-set 17 --out ", OUTPUTpath, "/pca_wild"))

system("./plink --bfile wild_standard_final --pca --chr-set 17 --out pca_wild") #run the pca code, specified the number of chromosomes


pca.data.wild <- read.table("./data/results/pca/pca_wild.eigenvec", header=FALSE) #load pca results from the eigenvec file
colnames(pca.data.wild)[2] <- "ID" #rename second column to allow for appending pop assignment later
eigenvalues <- read.table("./data/results/pca/pca_wild.eigenval", header=FALSE)$V1 #read in eigenvalues to calc variance
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