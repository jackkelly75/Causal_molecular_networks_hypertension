library(MRPC)
library(igraph)
library("AnnotationDbi") 
library("org.Hs.eg.db")
library(stringr)
set.seed(120)

#import the data
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
##import the adjusted data to find surrogate variables
load("covariate.Rdata") #mod1, covariate_data

#load adjusted expression data
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
load("eqtl_exprs_all_nobatch.Rdata")

#load eQTL data
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
load("eQTL_data.Rdata") #sig_eqtl, exprs
rm(exprs)

setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clustering/nobatch")
load("adjusted_module_info.Rdata") #moduleLabels, moduleColors, MEs
rm(MEs)
rm(moduleLabels)
colors <- unique(moduleColors)

#get list of all files for probes
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes/single_probes")
files <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
#import the colnames from the individual chr file
#columns <- read.table("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/chr2.txt", nrows = 1,comment.char = "")
columns <- read.table("/mnt/bmh01-rds/xiaoguang-cells/eQTL/data/chr2.txt", nrows = 1,comment.char = "")
columns <- as.character(columns)
# get the covariate data to attach to each table
condition <- as.character(covariate_data["Hypertension",])
names(condition) <- names(covariate_data["Hypertension",])
condition <- gsub("N" , 0, condition)
condition <- gsub("Y" , 1, condition)
hypertension <- as.numeric(condition)
names(hypertension) <- names(condition)
hypertension <- hypertension[colnames(eqtl_exprs)]
#get surrogate variables
#surrogate <- t(svd$sv)
#surrogate <- t(svd_leek$sv) #for leek adjustment
# leek adjustment may be better as larger number of surrogate variables may capture too much info
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0808-5
#rownames(surrogate) <- sprintf("surrogate%d", 1:nrow(surrogate))



for(i in colors){
	print(paste0("doing ", i))
	print(paste0("module ", which(colors == i), " of ", length(colors)))
	genes <- names(moduleColors[moduleColors == i])
	#from clumped file, keep the genotype info for the SNPs assocaited with each probe.
	#Go into each probe file and then from there get the SNPs
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes/single_probes")
	chr_files <- paste("./", genes, sep = "")
	chr_files <- intersect(chr_files, files)
	#length(chr_files)
	setwd(chr_files[1]) #set wd to the gene
	genotype <- read.table('clumped_individual_SNPs.csv', header = F, sep = "\t", dec = ".", comment.char = "", skip = 1)
	#relabel the columns as are incorrect in the file
	colnames(genotype) <- columns
	rownames(genotype) <- genotype[,3]
	genotype <- genotype[,10:ncol(genotype)]
	samples <- intersect(colnames(genotype),colnames(eqtl_exprs))
	genotype <- genotype[,samples]
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes/single_probes")
	for(x in 2:length(chr_files)){
		setwd(chr_files[x])
		indiv_geno <- read.table('clumped_individual_SNPs.csv', header = F, sep = "\t", dec = ".", comment.char = "", skip = 1)
		colnames(indiv_geno) <- columns
		rownames(indiv_geno) <- indiv_geno[,3]
		indiv_geno <- indiv_geno[!(rownames(indiv_geno) %in% rownames(genotype)),]
		indiv_geno <- indiv_geno[,10:ncol(indiv_geno)]
		samples <- intersect(colnames(indiv_geno),colnames(eqtl_exprs))
		indiv_geno <- indiv_geno[,samples]
		genotype <- rbind(genotype, indiv_geno)
		setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes/single_probes")
		print(paste0("completed ", x))
	}
	##############
	#Map to symbol
	##############
	exprs <- eqtl_exprs[genes,]
	genes_network <- rownames(exprs)
	symbols <- mapIds(org.Hs.eg.db, keys = genes_network, keytype = "ENSEMBL", column="SYMBOL")
	#if doesn't map, keep it as the ensembl rather than replace with NA
	for (p in 1:length(symbols)){
		if (is.na(symbols[p])) {
			symbols[p] <- names(symbols)[p]
		}
	}
	#duplicate genes are already removed
	rownames(exprs) <- symbols #relabel the genes to IDs

	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/nobatch")
	dir.create(i)
	setwd(i)
	save(exprs, genotype, file  = "genotype_exprs.Rdata")
	#ensure all in the same order
	#surrogate_variables <- surrogate[,colnames(genotype)]
	exprs <- exprs[,colnames(genotype)]
	hypertension <- hypertension[colnames(genotype)]
	#MRPC_data <- t(rbind(genotype, exprs, surrogate_variables, hypertension))
	MRPC_data <- t(rbind(genotype, exprs, hypertension))
	colnames(MRPC_data)[ncol(MRPC_data)] <- "HYPERTENSION"
	labels = colnames(MRPC_data)
	n <- nrow(MRPC_data)
	suffStat_C1 <- list(C = cor(MRPC_data), n = n)
	GV = sum(str_count(colnames(MRPC_data), pattern = "_b37"))
	save(MRPC_data, suffStat_C1, labels, file= "MRPC.Rdata")
	print(paste("completed ", i))
	#rm(surrogate_variables)
	rm(exprs)
	rm(genotype)
	rm(MRPC_data)
	rm(suffStat_C1)
	rm(labels)
}