#paths (edit according your directory set up)
export VCFTOOLS="/path/local/bin" #directory to your executable vcftools 
export PATH="/c/msys64/ucrt64/bin:$PATH" #to make it work in the console in R on a Windows computer; mileage may vary with MAC
mkdir -p "./results/diversity"

zcat ./data/processed/wild_standard_final.vcf.gz | "${VCFTOOLS}"/vcftools --vcf - --het --out ./results/diversity/wild_het

#command for producing windowed pi with vcftools, done for both populations:
zcat ./data/processed/refuge_standard_postsubset.vcf.gz | "${VCFTOOLS}"/vcftools --vcf - --window-pi 50000 --out ./results/diversity/refuge-nucleotide
zcat ./data/processed/ranch_standard_postsubset.vcf.gz | "${VCFTOOLS}"/vcftools --vcf - --window-pi 50000 --out ./results/diversity/ranch-nucleotide
