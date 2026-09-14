##admixture plotting

#using LD pruned dataset: wild_LDpruned_05

#setting working directory to folder with admixture
cd /Users/tylerbostwick/Documents/Masters_Work/Analyses/ADMIXTURE/dist/admixture_macosx-1.3.0
#running admixture
  for k in {2..5}; do ./admixture --cv=10 wild_LDpruned_05.bed $k > wild_LDpruned_05k${k}.txt; done

#runs the admixture program for k 2-5, generating cross-validation, and writing it to a txt file for each k ran