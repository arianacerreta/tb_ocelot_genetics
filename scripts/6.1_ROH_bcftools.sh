#used bcftools to select roh, uses a HMM to select runs
#code run in termainl is as follows:
#within your bcftools directory

export LDPRUNED_WILD="/c/Users/kualc095/OneDrive - Texas A&M University - Kingsville/Ocelots/LEPA_Manuscript1_Bostwick/tb_ocelot_genetics/data/processed/addtl_filter" #adjust path accordingly
#clean the header of the vcf for use in bcftools- removes a flag placed by plink that causes errors
./bcftools view -h "${LDPRUNED_WILD}"/roh_LDpruned_05.vcf.gz | grep -v "##chrSet" > "${LDPRUNED_WILD}"/clean_header.txt
./bcftools reheader -h "${LDPRUNED_WILD}"/clean_header.txt "${LDPRUNED_WILD}"/roh_LDpruned_05.vcf.gz -o "${LDPRUNED_WILD}"/clean_roh_LDpruned_05.vcf.gz

#skipping reheading
#make dir for roh results
export ROH_DIR="/c/Users/kualc095/OneDrive - Texas A&M University - Kingsville/Ocelots/LEPA_Manuscript1_Bostwick/tb_ocelot_genetics/results/roh"
mkdir "$ROH_DIR"

#identify roh:
./bcftools roh -G30 --estimate-AF - --rec-rate 1.1e-8 -Or -o "${ROH_DIR}"/roh_LD05.txt \
"${LDPRUNED_WILD}"/clean_roh_LDpruned_05.vcf.gz
#key changes, now estimates allele frequencies from data instead of using default, uses the domestic cat recombination rate
#and is using LD pruned data to better fit model assumptions
#as of 4/29/26 -- fixed the rec-rate flag
#kept G30, lowered it down to G20 and had minimal change in percentages, no break in long runs with the correct rec rate
#bcftools documentation states "safe value to use 30 to account for GT errors"
#--estimate-AF -   #this means recalculate INFO/AC and INFO/AN on fly with all samples ("-"); in case legacy AN and AC were carried over