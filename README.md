## *Genomic effects of rare gene flow between inbred populations of ocelots (Leopardus pardalis) in the United States*

#### Tyler Bostwick et al. (add upon submission)

##### Please contact the first author for questions about the code or data: Tyler Bostwick (add email)
##### Secondary contact: Lisanne Petracca (Lisanne.Petracca@tamuk.edu)

_______________________________________________________________________________________

## Abstract

(Add Abstract here) 

### Table of Contents 

(Add detail on what is in each folder within the repo - see the template for all possible folders, but at a minimum you should have the folders below, see examples below.) 

### [Scripts](./scripts)

Contains scripts to run all analyses. 
 
### [Data](./data) 

Contains raw and processed data.

### [Results](./results)

Contains raw and processed results.  

### [Figures](./figures)

Contains pdf versions of all figures in manuscript. 

### Required Packages and Versions Used 

(here, list all required packages and the version you used, see examples) 

here_1.0.1

dplyr_1.0.5

### Details of Article 

(Citation here) 

### How to Use this Repository 

#### Filtering steps

1.[1_bcftools_preprocessing.sh](./scripts/1_bcftools_preprocessing.sh)

- select autosomes, remove indels and non-biallelic snps, filter for genotype quality >9 

2.[2_PLINK_preanalysis.R](./scripts/2_PLINK_preanalysis.R)

- rename chromosomes, give unique IDs to SNPS, subset to wild individuals, filtering 
