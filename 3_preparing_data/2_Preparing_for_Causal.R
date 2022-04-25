#############
# Preparing the files so that they can be used with causal analysis.

#module load apps/gcc/R/4.1.0
#get list of sig probes and list of directories
setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
sig_eqtl <- read.table('sorted_eqtl_results_p0.01.csv', header = TRUE, sep = ",", dec = ".",comment.char = "")
sig_probes <- unique(sig_eqtl[,1])

dirs <- list.dirs()
#select the dirs that are probes (the first 2 are ./ and clump, and the last is GeneID (should be removed at earlier step but forgot to, adapt this if I do delete it))
dirs <- dirs[3:(length(dirs)-1)]

#see which don't have any data in them (will reduce number of probes)
null_probe <- vector()
for(x in 1:length(dirs)){
  setwd(dirs[x])
  info = file.info('clumped_individual_SNPs.csv')
  null_probe[x] <- info$size
  gene_name <- read.table('associated_SNPs.csv', header = F, sep = ",", dec = ".",comment.char = "")
  names(null_probe)[x] <- as.character(gene_name[1,1])
  print(paste0("completed ", x))
  setwd("/mnt/iusers01/bk01/m20349jk/scratch/single_probes")
}

null_probe <- names(null_probe[null_probe == 0])
#48 probes to remove
#null_probe_file <- paste("./", null_probe, sep = "")
#dirs <- dirs[!(dirs %in% null_probe_file)]
sig_probes <- sig_probes[!sig_probes %in% null_probe]


# get the expression data and only keep the probes of interest
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data")
exprs <- read.table('gene_expression.txt', header = TRUE, dec = ".",comment.char = "", row.names = 4)
probe_chr <- exprs[,1:3] #get the first 3 columns which give chr number (keep col 2+3 to maintain the table format)
exprs <- exprs[,4:ncol(exprs)] # remove the first 4 columns which give gene position
exprs <- exprs[row.names(exprs) %in% sig_probes, ] #remove the probes with nothing in folder
nrow(exprs) # 17400 probes

#keep the sig_eqtls of the genes in exprs
sig_eqtl <- sig_eqtl[sig_eqtl[,1] %in% row.names(exprs),]

setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
save(sig_eqtl, exprs,probe_chr, file = "eQTL_data.Rdata")
saveRDS(sig_eqtl, file = "sig_eqtl.rds")
saveRDS(exprs, file = "exprs.rds")
saveRDS(probe_chr, file = "probe_chr.rds")