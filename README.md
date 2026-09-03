## *Genomic effects of rare gene flow between inbred populations of ocelots (Leopardus pardalis) in the United States*

#### Tyler Bostwick et al. (add upon submission)

##### Please contact the first author for questions about the code or data: Tyler Bostwick (add email)
##### Secondary contact: Lisanne Petracca (Lisanne.Petracca@tamuk.edu)

_______________________________________________________________________________________

## Abstract

(Add Abstract here) 

### Table of Contents 

### [Scripts](./scripts)

Contains scripts to run all analyses. Analyses were primarily conducted in R, however some shell code is also included when necessary. 
 
### [Data](./data) 

Contains input files, which are small files used to subset data. Raw data will need to be downloaded from [Zenodo](./data) and stored locally on your device due to file size constraints. If you fork this repository and run the code, it will create a directories for processed data locally. These folders are not pushed to the GitHub repository due to file size. 

NOTE: Ariana fix Zenodo link once data files uploaded

### [Results](./results)

Contains raw and processed results.  

### [Figures](./figures)

Contains pdf versions of all figures in manuscript. 

### Required Packages and Programs and Versions Used 

R packages:

here_1.0.1 (update)

dplyr_1.0.5 (update)

Additional Programs:

[bcftools](https://samtools.github.io/bcftools/howtos/index.html) v1.24

[PLINK](https://www.cog-genomics.org/plink/1.9/) v1.9.0-b.7.11

[PLINK](https://www.cog-genomics.org/plink/2.0/) v2.0.0-a.7.4

### Details of Article 

Bostwick, T et al. Genomic effects of rare gene flow between inbred populations of ocelots (*Leopardus pardalis*) in the United States. In Prep  

### How to Use this Repository 

#### Filtering steps

[1_bcftools_preprocessing.sh](./scripts/1_bcftools_preprocessing.sh)

- select autosomes, remove indels and non-biallelic snps, filter for genotype quality >9 

[2_PLINK_preanalysis.R](./scripts/2_PLINK_preanalysis.R)

- rename chromosomes, give unique IDs to SNPS, subset to wild individuals, filtering

#### Analyses

[3_PCA_analysis.R](./scripts/3_PCA_analysis.R)

- perform PCA on wild ocelots and plot

[4_KING-robust_analysis.R](./scripts/4_KING-robust_analysis.R)

- calls PLINK2 to calculate KING-robust kinship estimator for WILD individuals and plots
