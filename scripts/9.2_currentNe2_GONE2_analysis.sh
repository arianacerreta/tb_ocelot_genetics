#If running in docker on a Mac:
##Set up an image called ne-tools with both currentne2, gone2, gone2_tui installed and compiled
##set up a folder called ne_analysis where results will be written

#currentNe2
#ROH and kinship filters; no MAF, miss 90, biallelic, coverage depth 7, genotype quality 9
#add --mac 1 to filter out polymorphic loci since GONE2 already autofilters them, but I think something weird is going on
###with structure population flag

docker run --rm -v ~/ne_analysis:/data ne-tools currentne2 -r 1.1 -x -o currentNe2_r1.1_x_wild_kin_roh_filt.txt wild_kin_roh_filter.vcf

###with no structure population flag
docker run --rm -v ~/ne_analysis:/data ne-tools currentne2 -r 1.1 -o currentNe2_r1.1_no-x_wild_kin_roh_filt.txt wild_kin_roh_filter.vcf

#GONE2
docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_update_r1_1_x_E_M0_wild_kin_roh_filter_mac1 -E -M 0.00 wild_kin_roh_filter_mac1.vcf
###Loading data files
###Reading the input file/s took 1.12 sec
###Measuring d²
###Resampled -x: 10 converged round(s), 99 usable bin(s).
###Fitting Ne(t)
###Combo kept: L2 (bin-weighted d² residual 0.000001)

#maf0.01
docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_update_r1_1_x_E__M01_wild_kin_roh_filter_mac1 -E -M 0.01 wild_kin_roh_filter_mac1.vcf
###Loading data files
###Reading the input file/s took 0.43 sec
###Measuring d²
###Resampled -x: 10 converged round(s), 99 usable bin(s).
###Fitting Ne(t)
### Combo kept: L2 (bin-weighted d² residual 0.000001)


docker run --rm -v ~/ne_analysis:/data ne-tools gone2 -g 0 -x -r 1.1 -o gone2_update_r1_1_x_E_M05_wild_kin_roh_filter_mac1 -E -M 0.05 wild_kin_roh_filter_mac1.vcf
###Loading data files
###Reading the input file/s took 0.45 sec
###Measuring d²
###Resampled -x: 10 converged round(s), 99 usable bin(s).
###Fitting Ne(t)
### Combo kept: L2 (bin-weighted d² residual 0.000001)
