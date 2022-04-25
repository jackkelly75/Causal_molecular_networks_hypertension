
# get the clumped SNPs and extract individual level data for these

#!/bin/bash --login
#$ -cwd               # Run job from current directory
cd /mnt/iusers01/bk01/m20349jk/scratch/single_probes/
#cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes
for OUTPUT in $(printf "%s\n" ENSG* )
do
  cd ${OUTPUT}
  awk -F, '{print $4}' clumped_associated_SNPs.csv  > clumped_SNPs.csv
  #awk -F, '{print $4}' clumped_associated_SNPs_0.01.csv  > clumped_SNPs_0.01.csv
  fgrep -w -f clumped_SNPs.csv individual_SNPs.csv > clumped_individual_SNPs.csv #get list of indiv SNPs data for the clumped SNPs
  sed -i "1s/^/$(head -n1 individual_SNPs.csv)\n/" clumped_individual_SNPs.csv #add heading to the SNP table
  #fgrep -w -f clumped_SNPs_0.01.csv individual_SNPs.csv > clumped_0.01_individual_SNPs.csv #get list of indiv SNPs data for the clumped SNPs
  #sed -i "1s/^/$(head -n1 individual_SNPs.csv)\n/" clumped_0.01_individual_SNPs.csv #add heading to the SNP table
  echo "${OUTPUT} completed"
  cd /mnt/iusers01/bk01/m20349jk/scratch/single_probes/
  #cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes
done
