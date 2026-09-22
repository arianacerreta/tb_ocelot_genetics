## Genomic effects of rare gene flow between inbred populations of ocelots (*Leopardus pardalis*) in the United States

#### Bostwick, T.A., Cerreta, A.L., DeYoung, R.W., Smith, M.M., Martin A.M., Reeves, A.R., and L.S. Petracca

##### Please contact the first author for questions about the code or data: Tyler Bostwick (add email)
##### Secondary contact: Lisanne Petracca (Lisanne.Petracca@tamuk.edu)

_______________________________________________________________________________________

## Abstract

(Add Abstract here) 

### Table of Contents 

### [Scripts](./scripts)

Contains scripts to run all analyses. Analyses were primarily conducted in R, however some shell code is also included when necessary. 
 
### [Data](./data) 

Contains input files, which are small files used to subset data. Raw data will need to be downloaded from [Zenodo](./data) and stored locally on your device due to file size constraints. If you fork this repository and run the code, it will create directories for processed data locally. These folders are not pushed to the GitHub repository due to file size. You will need to download the ocelot genome, [GCA_058742865.1](https://ncbi.nlm.nih.gov/datasets/genome/GCA_058742865.1/), for normalization.  

NOTE: Ariana fix Zenodo link once data files uploaded

### [Results](./results)

Contains raw and processed results.  

### [Figures](./figures)

Contains pdf versions of all figures in manuscript. 

### Required Packages, Programs, and Versions Used 

R packages:

dplyr_1.2.1

ggplot2_4.0.3


Additional Programs:

[bcftools](https://samtools.github.io/bcftools/howtos/index.html) v1.24

[PLINK](https://www.cog-genomics.org/plink/1.9/) v1.9.0-b.7.11

[PLINK](https://www.cog-genomics.org/plink/2.0/) v2.0.0-a.7.4

[ADMIXTURE](https://dalexander.github.io/admixture/) v1.3.0

### Details of Article 

Bostwick, T.A., Cerreta, A.L., DeYoung, R.W., Smith, M.M., Martin A.M., Reeves, A.R., and L.S. Petracca. Genomic effects of rare gene flow between inbred populations of ocelots (*Leopardus pardalis*) in the United States. In prep.

### How to Use this Repository 

Fork this repository to your computer. Download additional input files from [Zenodo](./data). Download [ocelot reference genome](https://ncbi.nlm.nih.gov/datasets/genome/GCA_058742865.1/). To recreate all intermediate files begin with [1_bcftools_preprocessing.sh](./scripts/1_bcftools_preprocessing.sh) and follow throught the pipeline as listed below. To avoid long processing times in pre-processing steps, begin with intermediate file ```allInd_SNPs_autosomes_bi_gq9.vcf.gz```and [2_PLINK_preanalysis.R](./scripts/2_PLINK_preanalysis.R). 

*Most* paths are internally referenced within the code with notable exceptions being your path to bcftools, PLINK, ADMIXTURE, and original input files downloaded from Zenodo. Paths that will need to be updated are annotated within the code.

Ariana: update Zenodo link once created

#### Filtering steps

[1_bcftools_preprocessing.sh](./scripts/1_bcftools_preprocessing.sh)

- select autosomes, remove indels and non-biallelic snps, filter for genotype quality >9 

[2_PLINK_preanalysis.R](./scripts/2_PLINK_preanalysis.R)

- rename chromosomes, give unique IDs to SNPs, subset to wild individuals, filtering

#### Analyses

[3_PCA_analysis.R](./scripts/3_PCA_analysis.R)

- perform PCA on wild ocelots and plot

[4_KING-robust_analysis.R](./scripts/4_KING-robust_analysis.R)

- calls PLINK2 to calculate KING-robust kinship estimator for WILD individuals and plots

[5.1_admixture_analysis.sh](./scripts/5.1_admixture_analysis.sh)

- calls and runs ADMIXTURE

[5.2_admixture_plotting.R](./scripts/5.2_admixture_plotting.R)

- reads in data from ADMIXTURE runs for plotting

[6.1_ROH_bcftools.sh](./scripts/6.1_ROH_bcftools.sh)

- calls bcftools to perform ROH analysis

[6.2_ROH_plotting.R](./scripts/6.2_ROH_plotting.R)

- plots results from ROH analysis
