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
save(mod1, covariate_data, file = "covariate.Rdata")



#Expression data
exprs_filt <- exprs[,rownames(mod1)]
data <-  sapply(exprs_filt, as.numeric)
rownames(data) <- rownames(exprs)
save(data, file = "data.Rdata")

#" The analysis was adjusted for age, body mass index, sample ischemic time, sample sequencing center,
#and a variable number of surrogate variables. The number of surrogate variables calculated per tissue
#was 10% of the number of samples tested. Surrogate variables were calculated by SVASeq23 on the middle
#80% of expressed genes in each tissue. The top and bottom 10% were excluded to reduce the impact of highly
#and weakly expressed genes,23 which are more likely to generate artifactual read count values."
#https://pubmed.ncbi.nlm.nih.gov/31644355/

#use combat rather than combat_seq as we do not need to maintain count data
##adjust for sex
sex = as.vector(as.character(covariate_data["SexM",]))
sex = as.factor(sex)
data_adj_sex <- ComBat(data, batch = sex, mod = mod1, par.prior=TRUE, prior.plots=FALSE)


#adjust for transplant status
transplant = as.vector(as.character(covariate_data["transplant",]))
transplant = as.factor(transplant)
data_adj_trans <- ComBat(data, batch = transplant, mod = mod1, par.prior=TRUE, prior.plots=FALSE)

age = as.vector(as.character(covariate_data["Age",]))
age = as.factor(age)
data_adj <- ComBat(data_adj_trans, batch = age, mod = mod1, par.prior=TRUE, prior.plots=FALSE)

save(data_adj, file = "data_adj.Rdata")



##using the limma approach
##adjust for sex
sex = as.vector(as.character(covariate_data["SexM",]))
sex = as.factor(sex)
transplant = as.vector(as.character(covariate_data["transplant",]))
transplant = as.factor(transplant)
age <- as.vector(as.numeric(covariate_data["Age",]))
data_adj <- removeBatchEffect(data, batch=sex, batch2=transplant, covariates=age, design=mod1)

save(data_adj, file = "data_adj_limma.Rdata")
