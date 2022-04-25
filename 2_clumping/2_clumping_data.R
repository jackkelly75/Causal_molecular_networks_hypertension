#module load apps/gcc/R/4.1.0
library(data.table)
library(dplyr)
library(TwoSampleMR)

#dirs <- list.dirs(path = "/mnt/iusers01/bk01/m20349jk/scratch/single_probes/", full.names = TRUE, recursive = TRUE)
dirs <- list.dirs(path = "/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes", full.names = TRUE, recursive = TRUE)


dirs <- dirs[2:length(dirs)]
#testing time to run 100
dirs <- dirs[1:100]
ptm <- proc.time()

vec <- numeric(0)
orig_vec <- numeric(0)
n = 1
for (p in 1:length(dirs)){
	setwd(dirs[p])
	x <- fread(file = "associated_SNPs.csv") #read in file
	orig_vec[n] <- nrow(x)
	names(orig_vec)[n] <- x[1,1]
	data <- x[,c(2, 3, 4, 5, 7, 8, 9)] #get the columns of interest for clumping
	colnames(data) <- c("SNP", "chr_name", "chrom_start", "chr_pos", "pval", "beta", "se")
	imported <- TwoSampleMR::format_data(data)
	clumped <- TwoSampleMR::clump_data(
	  imported,
	  clump_kb = 10000,
	  clump_r2 = 0.1,
	  clump_p1 = 1,
	  clump_p2 = 1,
	  pop = "EUR"
	)
	clumped_associated_SNPs <- data[data$SNP %in% clumped[,1],]
	write.table(clumped_associated_SNPs, file = "clumped_associated_SNPs.csv" , col.names = F, row.names = F, quote = F, sep = ",")
	vec[n] <- nrow(clumped_associated_SNPs)
	names(vec)[n] <- x[1,1]
	n = n + 1
}

#end timing
end <- proc.time() - ptm
end

mean(vec) #5.84
mean(orig_vec) #153.36




########################## bash code
#mkdir clump
#cd clump


#clump the data in batches of 1000 (should take less than 2 hours each)
#18 files to be created to do this
#gedit clumping_01_1.txt
#gedit clumping_01_2.txt
#gedit clumping_01_3.txt
#gedit clumping_01_4.txt
#gedit clumping_01_5.txt
#gedit clumping_01_6.txt
#gedit clumping_01_7.txt
#gedit clumping_01_8.txt
#gedit clumping_01_9.txt
#gedit clumping_01_10.txt
#gedit clumping_01_11.txt
#gedit clumping_01_12.txt
#gedit clumping_01_13.txt
#gedit clumping_01_14.txt
#gedit clumping_01_15.txt
#gedit clumping_01_16.txt
#gedit clumping_01_17.txt
#gedit clumping_01_18.txt




##############create this script
#!/bin/bash --login
#$ -cwd               # Run job from current directory
## We now recommend loading the modulefile in the jobscript. Change the version as needed.
module load apps/gcc/R/4.1.0
R CMD BATCH  clumping_01_1.R  clumping_01_1.R.o$JOB_ID





#gedit clumping_01_1.R
#gedit clumping_01_2.R
#gedit clumping_01_3.R
#gedit clumping_01_4.R
#gedit clumping_01_5.R
#gedit clumping_01_6.R
#gedit clumping_01_7.R
#gedit clumping_01_8.R
#gedit clumping_01_9.R
#gedit clumping_01_10.R
#gedit clumping_01_11.R
#gedit clumping_01_12.R
#gedit clumping_01_13.R
#gedit clumping_01_14.R
#gedit clumping_01_15.R
#gedit clumping_01_16.R
#gedit clumping_01_17.R
#gedit clumping_01_18.R
