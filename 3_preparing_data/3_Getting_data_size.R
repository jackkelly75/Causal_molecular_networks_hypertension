#module load apps/gcc/R/4.1.0
library("AnnotationDbi") 
library("org.Hs.eg.db")
library(stringr)

#eQTL data size
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
load("eQTL_data.Rdata")
sink('data_size_overview.txt')
setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
#get list of all folders in the file
folders <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)

# the header in the "clumped_individual_SNPs.csv" files are sep differently to the other rows.
#All rows are tab sep, however the header is weird combination of spaces
#Get new headers to use for the files
setwd("/net/scratch2/m20349jk/single_probes/ENSG00000164104") # select a random file
string <- read.table('clumped_individual_SNPs.csv', sep = "\t", header = F, nrows = 1) #import the first line of the file with weird spacing
string <- gsub("([0-9])([a-z])", "\\1 \\2", string) #split out any of the samples that have no whitespace between them. Adds space between 
string <- gsub("\\s+", " ", str_trim(string)) #replace any number of spaces with just one space
col_names <- strsplit(string, " +")[[1]] # break spaces into seperate strings in a vector



for(i in 1:22){
	#print which chromome working with
	print(paste0("Info for Chromosome ", i))
	#get list of probes for the chromome
	probe_chr21 <- rownames(probe_chr[probe_chr[,1] == paste0("chr", i),])
	exprs_chr21 <- exprs[rownames(exprs) %in% probe_chr21,]
	print(paste0("From all expression data, there are ", nrow(exprs_chr21), " probes"))
	#get the files that contain chr1 probes
	setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
	chr21_files <- paste("./", rownames(exprs_chr21), sep = "")
	chr21_files <- intersect(chr21_files, folders)
	print(paste0("Of the 17400 sig probes ", length(chr21_files), " are from chromosome ", i))
	#from clumped file, keep the genotype info for the SNPs assocaited with each probe.
	#Go into each probe file and then from there get the SNPs
	setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
	setwd(chr21_files[1])
	genotype <- read.table('clumped_individual_SNPs.csv', header = F, sep = "\t", dec = ".",comment.char = "", skip = 1)
	colnames(genotype) <- col_names #don't include the first header which is chr and has been imported as rownames
	samples <- intersect(colnames(genotype),colnames(exprs_chr21))
	genotype <- genotype[,samples]
	setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
	for(x in 2:length(chr21_files)){
		setwd(chr21_files[x])
		indiv_geno <- read.table('clumped_individual_SNPs.csv', header = F, sep = "\t", dec = ".",comment.char = "", skip = 1)
		colnames(indiv_geno) <- col_names
		#indiv_geno <- indiv_geno[!(rownames(indiv_geno) %in% rownames(genotype)),]
		indiv_geno <- indiv_geno[,10:ncol(indiv_geno)]
		samples <- intersect(colnames(indiv_geno),colnames(exprs_chr21))
		indiv_geno <- indiv_geno[,samples]
		genotype <- rbind(genotype, indiv_geno)
		setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
	}
	#get the expression data for the chr21 files we are interested in
	exprs_chr21 = exprs_chr21[substr(chr21_files, start = 3, stop = 50),]
	print(paste0("There are ", nrow(exprs_chr21), " genes from chromosome ", i))
	print(paste0("There are ", nrow(genotype), " sig. SNPs from chromosome ", i))
	#Map to symbol
	genes_network <- rownames(exprs_chr21)
	symbols <- mapIds(org.Hs.eg.db, keys = genes_network, keytype = "ENSEMBL", column="SYMBOL")
	for (p in 1:length(symbols)){
		if (is.na(symbols[p])) {
			symbols[p] <- names(symbols)[p]
		}
	}
	print(paste0("There are ", length(unique(symbols)), " unique Genes"))
}

sink()






[1] "Info for Chromosome 1"
[1] "From all expression data, there are 1725 probes"
[1] "Of the 17400 sig probes 1725 are from chromosome 1"
[1] "There are 1725 genes from chromosome 1"
[1] "There are 9776 sig. SNPs from chromosome 1"
[1] "There are 1723 unique Genes"
[1] "Info for Chromosome 2"
[1] "From all expression data, there are 1222 probes"
[1] "Of the 17400 sig probes 1222 are from chromosome 2"
[1] "There are 1222 genes from chromosome 2"
[1] "There are 7517 sig. SNPs from chromosome 2"
[1] "There are 1221 unique Genes"
[1] "Info for Chromosome 3"
[1] "From all expression data, there are 1029 probes"
[1] "Of the 17400 sig probes 1029 are from chromosome 3"
[1] "There are 1029 genes from chromosome 3"
[1] "There are 6273 sig. SNPs from chromosome 3"
[1] "There are 1029 unique Genes"
[1] "Info for Chromosome 4"
[1] "From all expression data, there are 689 probes"
[1] "Of the 17400 sig probes 689 are from chromosome 4"
[1] "There are 689 genes from chromosome 4"
[1] "There are 4501 sig. SNPs from chromosome 4"
[1] "There are 689 unique Genes"
[1] "Info for Chromosome 5"
[1] "From all expression data, there are 853 probes"
[1] "Of the 17400 sig probes 853 are from chromosome 5"
[1] "There are 853 genes from chromosome 5"
[1] "There are 5027 sig. SNPs from chromosome 5"
[1] "There are 852 unique Genes"
[1] "Info for Chromosome 6"
[1] "From all expression data, there are 940 probes"
[1] "Of the 17400 sig probes 940 are from chromosome 6"
[1] "There are 940 genes from chromosome 6"
[1] "There are 7341 sig. SNPs from chromosome 6"
[1] "There are 939 unique Genes"
[1] "Info for Chromosome 7"
[1] "From all expression data, there are 907 probes"
[1] "Of the 17400 sig probes 907 are from chromosome 7"
[1] "There are 907 genes from chromosome 7"
[1] "There are 6063 sig. SNPs from chromosome 7"
[1] "There are 905 unique Genes"
[1] "Info for Chromosome 8"
[1] "From all expression data, there are 630 probes"
[1] "Of the 17400 sig probes 630 are from chromosome 8"
[1] "There are 630 genes from chromosome 8"
[1] "There are 4026 sig. SNPs from chromosome 8"
[1] "There are 630 unique Genes"
[1] "Info for Chromosome 9"
[1] "From all expression data, there are 687 probes"
[1] "Of the 17400 sig probes 687 are from chromosome 9"
[1] "There are 687 genes from chromosome 9"
[1] "There are 4262 sig. SNPs from chromosome 9"
[1] "There are 687 unique Genes"
[1] "Info for Chromosome 10"
[1] "From all expression data, there are 688 probes"
[1] "Of the 17400 sig probes 688 are from chromosome 10"
[1] "There are 688 genes from chromosome 10"
[1] "There are 4455 sig. SNPs from chromosome 10"
[1] "There are 687 unique Genes"
[1] "Info for Chromosome 11"
[1] "From all expression data, there are 1011 probes"
[1] "Of the 17400 sig probes 1011 are from chromosome 11"
[1] "There are 1011 genes from chromosome 11"
[1] "There are 6059 sig. SNPs from chromosome 11"
[1] "There are 1011 unique Genes"
[1] "Info for Chromosome 12"
[1] "From all expression data, there are 954 probes"
[1] "Of the 17400 sig probes 954 are from chromosome 12"
[1] "There are 954 genes from chromosome 12"
[1] "There are 5580 sig. SNPs from chromosome 12"
[1] "There are 952 unique Genes"
[1] "Info for Chromosome 13"
[1] "From all expression data, there are 317 probes"
[1] "Of the 17400 sig probes 317 are from chromosome 13"
[1] "There are 317 genes from chromosome 13"
[1] "There are 2191 sig. SNPs from chromosome 13"
[1] "There are 316 unique Genes"
[1] "Info for Chromosome 14"
[1] "From all expression data, there are 609 probes"
[1] "Of the 17400 sig probes 609 are from chromosome 14"
[1] "There are 609 genes from chromosome 14"
[1] "There are 3872 sig. SNPs from chromosome 14"
[1] "There are 607 unique Genes"
[1] "Info for Chromosome 15"
[1] "From all expression data, there are 610 probes"
[1] "Of the 17400 sig probes 610 are from chromosome 15"
[1] "There are 610 genes from chromosome 15"
[1] "There are 3542 sig. SNPs from chromosome 15"
[1] "There are 608 unique Genes"
[1] "Info for Chromosome 16"
[1] "From all expression data, there are 860 probes"
[1] "Of the 17400 sig probes 860 are from chromosome 16"
[1] "There are 860 genes from chromosome 16"
[1] "There are 5408 sig. SNPs from chromosome 16"
[1] "There are 859 unique Genes"
[1] "Info for Chromosome 17"
[1] "From all expression data, there are 1032 probes"
[1] "Of the 17400 sig probes 1032 are from chromosome 17"
[1] "There are 1032 genes from chromosome 17"
[1] "There are 6126 sig. SNPs from chromosome 17"
[1] "There are 1031 unique Genes"
[1] "Info for Chromosome 18"
[1] "From all expression data, there are 291 probes"
[1] "Of the 17400 sig probes 291 are from chromosome 18"
[1] "There are 291 genes from chromosome 18"
[1] "There are 1983 sig. SNPs from chromosome 18"
[1] "There are 291 unique Genes"
[1] "Info for Chromosome 19"
[1] "From all expression data, there are 1234 probes"
[1] "Of the 17400 sig probes 1234 are from chromosome 19"
[1] "There are 1234 genes from chromosome 19"
[1] "There are 7735 sig. SNPs from chromosome 19"
[1] "There are 1234 unique Genes"
[1] "Info for Chromosome 20"
[1] "From all expression data, there are 458 probes"
[1] "Of the 17400 sig probes 458 are from chromosome 20"
[1] "There are 458 genes from chromosome 20"
[1] "There are 3188 sig. SNPs from chromosome 20"
[1] "There are 458 unique Genes"
[1] "Info for Chromosome 21"
[1] "From all expression data, there are 195 probes"
[1] "Of the 17400 sig probes 195 are from chromosome 21"
[1] "There are 195 genes from chromosome 21"
[1] "There are 1598 sig. SNPs from chromosome 21"
[1] "There are 195 unique Genes"
[1] "Info for Chromosome 22"
[1] "From all expression data, there are 459 probes"
[1] "Of the 17400 sig probes 459 are from chromosome 22"
[1] "There are 459 genes from chromosome 22"
[1] "There are 2890 sig. SNPs from chromosome 22"
[1] "There are 459 unique Genes"

