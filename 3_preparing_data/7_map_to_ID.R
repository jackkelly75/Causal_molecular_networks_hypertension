#module load apps/gcc/R/4.1.0
#module load tools/env/proxy
library(pcalg) 
set.seed(120)
library("AnnotationDbi") 
library("org.Hs.eg.db")


setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
#load("data_adj_sv.Rdata")
#load("data_adj_leek.Rdata")
#load("data_adj_custom.Rdata")
#load("data_adj_all.Rdata")
load("data_adj_all_nobatch.Rdata")

data_adj <- data_adj_all


#Map to symbol
genes_network <- rownames(data_adj)
symbols <- mapIds(org.Hs.eg.db, keys = genes_network, keytype = "ENSEMBL", column="SYMBOL")
#if doesn't map, keep it as the ensembl rather than replace with NA
for (p in 1:length(symbols)){
	if (is.na(symbols[p])) {
		symbols[p] <- names(symbols)[p]
	}
}
#any dublicate gene names, only the gene with the highest MAD is kept
MAD <- vector(mode="numeric", length=0)
for( i in 1:nrow(data_adj)){                
        MAD[i] <- mad(data_adj[i,1:ncol(data_adj)])
}
ExprsMAD <- data.frame(symbols, MAD, data_adj)
#put the table in order by gene name and then MAD
ExprsMAD = ExprsMAD[order(ExprsMAD[,1], abs(ExprsMAD[,2]), decreasing = TRUE), ]
#remove dupicates (automatically removes anything after first instance, since is ordered by MAD only keeps highest MAD)
ExprsUniq <- ExprsMAD[!duplicated(ExprsMAD$symbol),]
eqtl_exprs <- ExprsUniq[,3:ncol(ExprsUniq)]
#this leaves 17382 genes (removes 18 genes)

save(eqtl_exprs, file = "eqtl_exprs_all.Rdata")
save(eqtl_exprs, file = "eqtl_exprs_all_nobatch.Rdata")