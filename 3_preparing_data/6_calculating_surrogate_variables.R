#PEER isn't working for install anymore. Instead using SVA.
#info on SVA- https://academic.oup.com/nar/article/42/21/e161/2903156?login=true
#module load apps/gcc/R/4.1.0
#module load tools/env/proxy

#import packages
library(sva)
library(limma)
#import the data
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
##import the adjusted data to find surrogate variables
#load("data_adj.Rdata")
load("data_adj_limma.Rdata")
load("covariate.Rdata") #mod1, covariate_data

#inverse function of data as is raised log2 and don't want that for the svaseq - use sva instead
#data <- 2^data_adj


#create surrogate variables
#svd <- svaseq(data_adj,mod1)
#svd <- sva(data_adj,mod1)
#32 surrogate variables for combat
#31 surrogate variables for limma
#n.sv = num.sv(data_adj,mod1,method="leek", seed = 2)
#svd_leek <- sva(data_adj,mod1, n.sv=n.sv)
#identifies 2 SVs

#use prop_var_testing.R to identify the number of SVs to use
svd_custom <- sva(data_adj,mod1, n.sv=6)
#dat - 	The transformed data matrix with the variables in rows and samples in columns
#mod - The model matrix being used to fit the data


#rownames(svd$sv) <- colnames(data_adj)
#rownames(svd_leek$sv) <- colnames(data_adj)
rownames(svd_custom$sv) <- colnames(data_adj)


setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
#save(svd, svd_leek, file = "svd.Rdata")
save(svd_custom, file = "svd_custom.Rdata")

#data_adj_sv <- removeBatchEffect(data_adj, covariates=svd$sv, design=mod1)
#data_adj_leek <- removeBatchEffect(data_adj, covariates=svd_leek$sv, design=mod1)
data_adj_custom <- removeBatchEffect(data_adj, covariates=svd_custom$sv, design=mod1)

#save(data_adj_sv, file = "data_adj_sv.Rdata")
#save(data_adj_leek, file = "data_adj_leek.Rdata")
save(data_adj_custom, file = "data_adj_custom.Rdata")