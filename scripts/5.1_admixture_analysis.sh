##admixture plotting

#using LD pruned dataset: wild_LDpruned_05

#setting paths to folder with admixture
export ADMIXTURE="/path/admixture_macosx-1.3.0" #update directory to your executable admixture
export INPUT_FILE="./data/processed/addtl_filter/wild_LDpruned_05.bed" #no need to change if forked from GitHub
export OUTPUT_DIR="./results/admixture" #no need to change if forked from GitHub

#running admixture

for run in {1..2}; do
 mkdir -p "${OUTPUT_DIR}/run${run}"
 (
  cd "${OUTPUT_DIR}/run${run}" || exit 1
  for k in {1..5}; do 
  "${ADMIXTURE}"/admixture -s $RANDOM --cv=10 "$INPUT_FILE" $k > "${OUTPUT_DIR}/run${run}"/wild_LDpruned_05_k${k}_run${run}.txt; 
  done
)
done
#runs the admixture program for k 2-5, generating cross-validation, and writing it to a txt file for each k ran

#save CV errors for plotting
for run in {1..10}; do
  (
    cd "${OUTPUT_DIR}/run${run}" || exit 1
    grep "CV" *txt | awk -v run="${run}" '{print $3,$4, run}' | sed -e 's/(//;s/)//;s/://;s/K=//' > CV_per_K_run${run}.txt
  )
done
