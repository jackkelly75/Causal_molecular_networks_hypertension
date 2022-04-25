#get list of the modules that have connections to hypertension
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/nobatch")
col <- list.dirs(recursive = F)
modules <- vector()
number <- vector()
for(i in col){
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/nobatch")
	setwd(i)
	folder <- sub('./', '', i)
	print(paste0("Start ", folder))
	if (file.exists(paste0("MRPC_ancestor_adjacency_", folder, ".csv"))){
		adj <- read.csv(paste0("MRPC_ancestor_adjacency_", folder, ".csv"))
		if(dim(adj)[1] > 1){
			modules <- c(modules, i)
			number <- c(number, nrow(adj))
		}
	}
}

