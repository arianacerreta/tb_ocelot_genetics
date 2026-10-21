if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")

library(dartR)
library(vcfR)
NEestimator<-"F:/2_TB_Working_Files/NeEstimator"

vcf<-gl.read.vcf("./data/processed/ranch_standard_postsubset.vcf")

dir.create("./results/current_Ne")
ranch_ne<-gl.LDNe(vcf,
                  outfile = "ranchLDNe.txt",
                  outpath = "./results/current_Ne/",
                  neest.path = NEestimator,
                  critical = c(0.05),
                  singleton.rm = TRUE,
                  mating = "random",
                  verbose = 5)
