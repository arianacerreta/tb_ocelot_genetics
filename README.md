## Genomic effects of rare gene flow between inbred populations of ocelots (*Leopardus pardalis*) in the United States

#### T. A. Bostwick, A. L. Cerreta, R. W. DeYoung, M. M. Smith, L. Caua, B. Davis, J. Janecka, A. M. Martin, A. R. Reeves, M. E. Tewes, and L. S. Petracca. 

##### Please contact the first author for questions about the code or data: Tyler Bostwick (add email)
##### Secondary contact: Lisanne Petracca (Lisanne.Petracca@tamuk.edu)

_______________________________________________________________________________________

## Abstract

(Add Abstract here) 

### Table of Contents 

### [Scripts](./scripts)

Contains scripts to run all analyses. Analyses were primarily conducted in R, 
however some shell code is also included when necessary. 
 
### [Data](./data) 

Contains input files, which are small files used to subset data. Raw data will 
need to be downloaded from [Zenodo](./data) and stored locally on your device 
due to file size constraints. If you fork this repository and run the code, it 
will create directories for processed data locally. These folders are not pushed 
to the GitHub repository due to file size. You will need to download the ocelot genome,
[GCA_058742865.1](https://ncbi.nlm.nih.gov/datasets/genome/GCA_058742865.1/), for normalization.  

NOTE: Ariana fix Zenodo link once data files uploaded

### [Results](./results)

Contains most results. Admixture results are not pushed to the GitHub repo. They 
should be generated to your local directories using the scripts.

### [Figures](./figures)

Contains pdf versions of all figures in manuscript or supplementals. 

### Required Packages, Programs, and Versions Used 

R packages:

dplyr_1.2.1

ggplot2_4.0.3

tidyr_1.3.2

stringr_1.6.0

forcats_1.0.1

cowplot_1.2.0

Additional Programs:

[bcftools](https://samtools.github.io/bcftools/howtos/index.html) v1.24

[PLINK](https://www.cog-genomics.org/plink/1.9/) v1.9.0-b.7.11

[PLINK](https://www.cog-genomics.org/plink/2.0/) v2.0.0-a.7.4

[ADMIXTURE](https://dalexander.github.io/admixture/) v1.3.0

[vcftools](https://vcftools.github.io/man_latest.html) v0.1.17

### Details of Article 

Bostwick, T. A., A. L. Cerreta, R. W. DeYoung, M. M. Smith, L. Caua, B. Davis, 
J. Janecka, A. M. Martin, A. R. Reeves, M. E. Tewes, and L. S. Petracca. 
Genomic effects of rare gene flow between inbred populations of ocelots (*Leopardus pardalis*) in 
the United States. In prep.

### How to Use this Repository 

Fork this repository to your computer. Download additional input files from 
[Zenodo](./data). Download [ocelot reference genome](https://ncbi.nlm.nih.gov/datasets/genome/GCA_058742865.1/). 
To recreate all intermediate files begin with [1_bcftools_preprocessing.sh](./scripts/1_bcftools_preprocessing.sh) 
and follow through the pipeline as listed below. To avoid long processing times 
in pre-processing steps, begin with intermediate file ```allInd_SNPs_autosomes_bi_gq9.vcf.gz```and [2_PLINK_preanalysis.R](./scripts/2_PLINK_preanalysis.R). 

*Most* paths are internally referenced within the code with notable exceptions 
being your path to bcftools, PLINK, ADMIXTURE, vcftools, and original input files downloaded 
from Zenodo. Paths that will need to be updated are annotated within the code or are indicated with ```/path/```.

Ariana: update Zenodo link once created

#### Filtering steps

[1_bcftools_preprocessing.sh](./scripts/1_bcftools_preprocessing.sh)

- select autosomes, remove indels and non-biallelic snps, filter for genotype quality >9 

[2_PLINK_preanalysis.R](./scripts/2_PLINK_preanalysis.R)

- rename chromosomes, give unique IDs to SNPs, subset to wild individuals, filtering

#### Analyses

[3.1_bcftools_het.sh](./scripts/3.1_bcftools_het.sh)

- calculates individual heterozygosity, π, F with bcftools

[3.2_diversity_stats_plotting.R](./scripts/3.2_diversity_stats_plotting.R)

- uses outputs from ```3.1_bcftools_het.sh``` to summarize and plot data

[4_PCA_analysis.R](./scripts/4_PCA_analysis.R)

- perform PCA on wild ocelots and plot

[5_KING-robust_analysis.R](./scripts/5_KING-robust_analysis.R)

- calls PLINK2 to calculate KING-robust kinship estimator for WILD individuals and plots

[6.1_admixture_analysis.sh](./scripts/6.1_admixture_analysis.sh)

- calls and runs ADMIXTURE

[6.2_admixture_plotting.R](./scripts/6.2_admixture_plotting.R)

- reads in data from ADMIXTURE runs for plotting

[7.1_ROH_bcftools.sh](./scripts/7.1_ROH_bcftools.sh)

- calls bcftools to perform ROH analysis

[7.2_ROH_plotting.R](./scripts/7.2_ROH_plotting.R)

- plots results from ROH analysis
