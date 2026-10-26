#If running in docker on a Mac:
##Set up an image called ne-tools with both currentne2, gone2, gone2_tui installed and compiled
##set up a folder called ne_analysis where results will be written

###with structure population flag

docker run --rm -v ~/ne_analysis:/data ne-tools currentne2 -r 1.1 -x -o currentNe2_r1.1_x_wild_no_miss_nomaf.txt wild_Ne_no_miss_nomaf.vcf

###with no structure population flag
docker run --rm -v ~/ne_analysis:/data ne-tools currentne2 -r 1.1 -o currentNe2_r1.1_no-x_wild_no_miss_nomaf.txt wild_Ne_no_miss_nomaf.vcf

###check results—estimates from default (no -x flag) were <<< that with -x indicating enough pop structuring to warrant -x flag; 
###use -x flag results and use -x flag for GONE2

#GONE2; no maf
docker run --rm -it -v ~/ne_analysis:/data ne-tools gone2_tui -g 0 -x -r 1.1 -o gone2_r1_1_x_E_wild_no_miss_nomaf -E wild_Ne_no_miss_nomaf.vcf

#maf 0.01
docker run --rm -it -v ~/ne_analysis:/data ne-tools gone2_tui -g 0 -x -r 1.1 -o gone2_r1_1_x_E_MAF01_wild_no_miss_nomaf -E -M 0.01 wild_Ne_no_miss_nomaf.vcf

#maf 0.05
docker run --rm -it -v ~/ne_analysis:/data ne-tools gone2_tui -g 0 -x -r 1.1 -o gone2_r1_1_x_E_MAF05_wild_no_miss_nomaf -E -M 0.05 wild_Ne_no_miss_nomaf.vcf


#ROH and kinship filters; no MAF, miss 90, biallelic, coverage depth 7, genotype quality 9
#add --mac 1 to filter out polymorphic loci since GONE2 already autofilters them, but I think something weird is going on
###with structure population flag

docker run --rm -v ~/ne_analysis:/data ne-tools currentne2 -r 1.1 -x -o currentNe2_r1.1_x_wild_kin_roh_filt.txt wild_kin_roh_filter.vcf

###with no structure population flag
docker run --rm -v ~/ne_analysis:/data ne-tools currentne2 -r 1.1 -o currentNe2_r1.1_no-x_wild_kin_roh_filt.txt wild_kin_roh_filter.vcf

#GONE2
docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_r1_1_x_E_M0_wild_kin_roh_filter_mac1 -E -M 0.00 wild_kin_roh_filter_mac1.vcf

###Loading data files
###Reading the input file/s took 1.04 sec
###Measuring d²
###Resampled -x: 10 converged round(s), 99 usable bin(s).
###Fitting Ne(t)
###Combo kept: L2 (bin-weighted d² residual 0.000001)

#force it to use all SNPs
docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_r1_1_x_E_M0_wild_kin_roh_filter_mac1_allpoly -E -M 0.00 -s 477244 wild_kin_roh_filter_mac1.vcf
###this forces a single pass which really introduced variability in the estimates...
###Loading data files
###Reading the input file/s took 1.06 sec
###Measuring d²
###Single-pass -x (all SNPs and individuals analysed): 99 usable bin(s).
###Fitting Ne(t)
###Combo kept: L2 (bin-weighted d² residual 0.000001)

docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_r1_1_x_E__M01_wild_kin_roh_filter_mac1 -E -M 0.01 wild_kin_roh_filter_mac1.vcf
###—-ERROR: setting progress bigger than max value
###*but still runs...and produces an output similar to the others
###Final cmd output:
###Resampled -x: 10 converged round(s), 99 usable bin(s).
###Fitting Ne(t)
###Combo kept: L2 (bin-weighted d² residual 0.000001)
###---for some reason it did not default to 20000 when applying the MAF filter

docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_r1_1_x_E__M01_wild_kin_roh_filter_mac1_force20000 -E -M 0.01 -s 20000 wild_kin_roh_filter_mac1.vcf
###—-ERROR: setting progress bigger than max value
###Final cmd output:
###Resampled -x: 10 converged round(s), 99 usable bin(s).
###Fitting Ne(t)
###Combo kept: L2 (bin-weighted d² residual 0.000001)
###***this did not force it to use 20000 SNPs...it used the same as the run without the -s 20000 flag

docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_r1_1_x_E__M05_wild_kin_roh_filter_mac1 -E -M 0.05 wild_kin_roh_filter_mac1.vcf
###—-ERROR: setting progress bigger than max value
###Final cmd output:
###Resampled -x: 10 converged round(s), 99 usable bin(s).
###Fitting Ne(t)
###Combo kept: L2 (bin-weighted d² residual 0.000001)