#module load apps/gcc/R/4.1.0

#load packages
library(stringr)
#WD with the Rdata files that include the number of SNPs associated with each probe
setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
#setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes")

x <- list.files()
x <- str_subset(x, "per_probe", negate = FALSE)
orig_files <- str_subset(x, "orig", negate = FALSE)
vec_files <- str_subset(x, "orig", negate = TRUE)


all_orig_vec <- vector()
for(i in orig_files){
	load(i)
	all_orig_vec  <- c(all_orig_vec, orig_vec)
	rm(orig_vec)
}


all_vec <- vector()
for(i in vec_files){
	load(i)
	all_vec  <- c(all_vec, vec)
	rm(vec)
}


mean(all_orig_vec) #182.1133
mean(all_vec) #6.331167


