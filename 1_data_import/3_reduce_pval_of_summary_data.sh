
cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data

#reduce the pval cut off for eQTLs to be 1e-7
awk -F"," '(NR>1) && ($4 <= 1e-7) ' /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/summary_p0.01.csv > ./summary_p1e-7.csv

awk -F"," '(NR>1) && ($4 <= 1e-4) ' /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/summary_p0.01.csv > ./summary_p1e-4.csv

awk -F"," '(NR>1) && ($4 <= 1e-2) ' /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/summary_p0.01.csv > ./summary_p1e-2.csv
