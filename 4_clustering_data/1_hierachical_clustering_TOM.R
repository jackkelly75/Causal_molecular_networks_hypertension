library(pcalg) 
set.seed(120)
library("AnnotationDbi") 
library("org.Hs.eg.db")
library(nortest) #normal distribution test
library(WGCNA)
enableWGCNAThreads()

setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
load("eqtl_exprs.Rdata")

setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clustering/nobatch")
##########
#	TOM  #
##########
# Choose a set of soft-thresholding powers
powers = c(c(1:15), seq(from = 15, to=20, by=2))
sft = pickSoftThreshold(t(eqtl_exprs), powerVector = powers, verbose = 5, networkType = "signed")
#use threshold of 0.85
# Plot the results
#pdf("TOM_softThreshold.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of 0.85
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#set the softPower found best above
softPower = 12
##create adjacency matrix
adjacency = adjacency(t(eqtl_exprs), power = softPower, type = "signed")
save(adjacency, file = "adjacency.Rdata")
#######plot to see the scale free topology
k=softConnectivity(t(eqtl_exprs),power=softPower)
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
########create a topological overlap matrix
TOM = TOMsimilarity(adjacency, TOMType="signed")
rownames(TOM) <- rownames(adjacency)
colnames(TOM) <- colnames(adjacency)
save(TOM, file = "TOM.Rdata")
##########turn into distance matrix
dissTOM = 1-TOM
save(dissTOM, file = "dissTOM.Rdata")

