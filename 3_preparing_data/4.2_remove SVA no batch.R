#adjust for SVA no batch
#module load tools/env/proxy
#module load apps/gcc/R/4.1.0

#import packages
library(sva)
library(limma)

#import the data
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
load("eQTL_data.Rdata")
# sig_eqtl summary data
# exprs is expression data

#import the covariate imformation
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/")
x <- read.table("covariates.txt", header = T, row.names = 1)
y <- t(read.table("hypertension_phenotype_eQTL.csv", header = T, row.names = 1, sep = ","))
covariate_data <- rbind(x, y)
#remove any samples that are unknown hypertension
covariate_dis <- covariate_data[,covariate_data["Hypertension",] == "Y"]
covariate_con <- covariate_data[,covariate_data["Hypertension",] == "N"]
covariate_data <- cbind(covariate_dis, covariate_con)
#extract the hypertension data and make a model matrix of it
group = as.vector(as.character(covariate_data["Hypertension",]))
group = as.factor(group)
mod1 = model.matrix(~group)
rownames(mod1) <- colnames(covariate_data) # set row names to the samples

setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
#save(mod1, covariate_data, file = "covariate.Rdata")
#Expression data
exprs_filt <- exprs[,rownames(mod1)]
data <-  sapply(exprs_filt, as.numeric)
rownames(data) <- rownames(exprs)
#save(data, file = "data.Rdata")


#use combat rather than combat_seq as we do not need to maintain count data
#use prop_var_testing.R to identify the number of SVs to use
svd_custom <- sva(data,mod1, n.sv=5)
rownames(svd_custom$sv) <- colnames(data)

data_adj_all <- removeBatchEffect(data, covariates=svd_custom$sv, design=mod1)


setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
save(svd_custom, file = "svd_all_nobatch.Rdata")
save(data_adj_all, file = "data_adj_all_nobatch.Rdata")