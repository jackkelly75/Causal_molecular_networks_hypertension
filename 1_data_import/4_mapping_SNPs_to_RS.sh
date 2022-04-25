#Mapping the SNPs
#Use the UKBB GWAS imputed v3 (https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=178908679)
#The file with rsIDs and positions is downloaded from https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz

cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/map_to_rsId
#download the variants.tsv.bgz to this folder
gunzip -c variants.tsv.bgz > variants.tsv
#reduce number of columns for better viewing and faster loading
awk -v s=1 '{print $1,$2,$3,$4,$5,$6,$7 >"filtered_variants.tsv"}' variants.tsv

#get the SNPs that are present at 0.01 pvalue level
awk -F"," '{print $2}' /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/summary_p0.01.csv > SNP_col.csv
sed 's/....$//' < SNP_col.csv > SNP_col1.csv
sed 's/_/:/' SNP_col1.csv > SNP_info.csv
paste SNP_info.csv SNP_col.csv > eQTL_SNPs_combined.csv

#load R to run this below - very easy to run and good to look at the data to see if better than previous approach
module load apps/binapps/rstudio/1.1.463   
qrsh -l short -cwd -V rstudio
##########
##########
##########
library(data.table)
library(dplyr)
x <- fread(file = "filtered_variants.tsv", sep = " ") #list of variants from UKBB GWAS imputed v3
x <- x[,6:7] #extract only the RSid and matching variant IDs
x <- x[x$rsid %like% "rs", ] #only keep those that have rs ids
y <- fread(file = "eQTL_SNPs_combined.csv", sep = "\t", header= F) #import the SNPs at 0.01 pvalue level
colnames(y)[1] <- "varid"
y <- unique(y)
x <- unique(x)
l <- merge(x, y, by="varid") #match the SNPs and the rsIDs
print(paste0("original has ",nrow(y)))
#"original has 1961319"
print(paste0("merged has ", nrow(l)))
#"merged has 1739651"
l <- l[,2:3]
#export table with mapping between SNP IDs and rsIDs
write.table(l, file = "mapped_eQTLs.csv", sep = ",", row.names =F, col.names=F, quote = F)
##########
##########
##########

wc -l /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/summary_p0.01.csv
#3607115
#join the file with rsId and positions to the summary data file
join -t, -1 2 -2 2 <(sort -t ',' -k 2 mapped_eQTLs.csv) <(sort -t ',' -k 2 /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/summary_p0.01.csv) > rs_eqtl_results_p0.01.csv
wc -l rs_eqtl_results_p0.01.csv #summary data mapped to RSid
#down from 3607115 to 3172616 SNPs in summary data (SNPs are in multiple times as mapped to individual genes)


cut -d, -f1 rs_eqtl_results_p0.01.csv > temp.csv
sed 's/_/,/g' temp.csv > temp1.csv
paste -d, temp1.csv rs_eqtl_results_p0.01.csv > temp2.csv
awk -F, '{print $8","$7","$1","$2","$6","$9","$10","$11","$12}' temp2.csv > temp3.csv

( echo -e "GeneID,SNP,chr_name,chrom_start,SNP_ID,TSS_dist,pval,beta,se"; cat temp3.csv ) >  sorted_eqtl_results_p0.01.csv
rm temp*.csv


#file now available at 
#/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/map_to_rsId/sorted_eqtl_results_p0.01.csv
