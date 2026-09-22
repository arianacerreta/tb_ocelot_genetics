##admixture plotting

#using LD pruned dataset: wild_LDpruned_05

#setting paths to folder with admixture
export ADMIXTURE="/path/admixture_macosx-1.3.0" #update directory to your executable admixture
export INPUT_FILE="./data/processed/addtl_filter/wild_LDpruned_05.bed" #no need to change if forked from GitHub
export OUTPUT_DIR="./results/admixture" #no need to change if forked from GitHub

#make new output dir
mkdir -p "$OUTPUT_DIR"

#running admixture
  for k in {2..5}; do "${ADMIXTURE}"/admixture --cv=10 "$INPUT_FILE" $k > "${OUTPUT_DIR}"/wild_LDpruned_05k${k}.txt; done

#runs the admixture program for k 2-5, generating cross-validation, and writing it to a txt file for each k ran