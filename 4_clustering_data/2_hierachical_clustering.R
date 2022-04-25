library(pcalg) 
set.seed(120)
library(flashClust)
library(km2gcn)
library(WGCNA)
#enableWGCNAThreads()


setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/")
load("eqtl_exprs.Rdata")
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clustering/nobatch")
load("dissTOM.Rdata")


#######
#WGCNA hierarchical
#######
geneTree = flashClust(as.dist(dissTOM), method="average") #flashClust is faster hclust
minModuleSize = 20 #set the smallest module size possible
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=4, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
#deepSplit can be changed to make more or less clusters
dynamicColors= labels2colors(dynamicMods)
#merge any close modules
MEDissThres = 0.2
merge = mergeCloseModules(t(eqtl_exprs), dynamicColors, cutHeight = MEDissThres, verbose = 3)

# merged module colors
moduleColors = merge$colors
names(moduleColors) <- colnames(t(eqtl_exprs))
#eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#starting with grey make list of colours and assign numerical value to each color
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
names(moduleLabels) <- colnames(t(eqtl_exprs))
MEs=mergedMEs
rownames(MEs)<-rownames(t(eqtl_exprs))

################################
#print the modules out to csv  #
################################
#setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clustering/hier")
#colors <- unique(moduleColors)
#for (i in 1:length(colors)){
#	MyData <- names(which(moduleColors == colors[i]))
#	write.csv(MyData, file = paste(colors[i],"_TOM_hierachical.csv", sep = ""))
#}


##############
#k-means adjustment of clustering
rm(dissTOM)
rm(merge)
rm(geneTree)
load("TOM.Rdata")

net = list(moduleColors, MEs)
names(net) <- c("moduleColors", "MEs")
expr.data = t(eqtl_exprs)
rm(eqtl_exprs)
n.iterations=400
meg = 0
tom.matrix=TOM
plot.evolution=TRUE
plot.go=FALSE
debug=NULL
net.type="signed"
min.genes.for.grey=20

#Step 1. Let D be the expression data in which dij in D represents the expression value for
#sample i and gene j, being s samples and g genes in total
cat("Working with",nrow(expr.data),"samples and",ncol(expr.data),"genes\n")
##Working with 455 samples and 17382 genes

#Step 2. Construct the partition by the WGCNA process, let P_D={m_1, m_2, ..., m_n} be that partition where m_k is the k-th module.
cat("The network includes",length(net$moduleColors),"genes and ", length(unique(net$moduleColors)),"modules\n")
#The network includes 17382 genes and  152 modules (110 for 0.2 merge)
partition.in.colors <- net$moduleColors

#Step 3. Get the eigengenes for each module within the partition, E={e_1, e_2, ..., e_n}
if(sum(partition.in.colors == "grey") < min.genes.for.grey){
	eigengenes = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=TRUE)
}else{
	eigengenes = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=F)
}
cat("We got",length(eigengenes$eigengenes)," eigengene vectors\n")
#We got 21  eigengene vectors
centroid.labels <- substring(names(eigengenes$eigengenes),3)
print("Module colors are")
print(centroid.labels)
#Step 4. Set up the k-means clustering
	#Step 4.1. Set k to n
k <- length(eigengenes$eigengenes)
	#Step 4.2. Set the centroids C to the eigengenes E, thus C to E
createCentroidMatrix <- function(eigengenes){
  my.matrix <- NULL
  for(eigengene in eigengenes){
    my.matrix <- cbind(my.matrix,eigengene)
  }
  return(my.matrix)
}
centroids <- createCentroidMatrix(eigengenes$eigengenes)
colnames(centroids) <- colnames(eigengenes$eigengenes)

#Step 5. Run the algorithm and monitor its evolution
#Step 5.1 Set iterations to 0
#Step 5.2 Create a new partition P', given C with n modules such that, for each gene, 1 <= j <= g, g_j belongs to the module c_t in C such that a distance meassure d(g_j,c_t) is minimum.
#Step 5.3 Calculate eigengenes of P', giving a new E'
#Step 5.4 Evaluate the progress. If progress done, set iterations to iterations + 1 and C to E' and go to step 5.2
#Step 5.5 Finish
partitions <- list()
new.partition <- match(partition.in.colors, centroid.labels)
names(new.partition) <- centroid.labels[new.partition]
partitions[[1]] <- new.partition

#Launch the iterations
exchanged.genes = meg + 1
iteration = 1
getBestModuleCor <- function(gene,centroids,signed=TRUE){
  return(which.max(0.5 * (1 + cor(centroids, gene, use = "all.obs", method = "pearson"))))
}
getExchangedGenes <- function(old.partition,new.partition){   		stopifnot(length(old.partition) == length(new.partition))
	return(old.partition[old.partition != new.partition])
}
getNewCentroids <- function(expr.data,partition.in.colors,centroid.labels,mgg){
  if(sum(partition.in.colors == "grey") < mgg)
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=TRUE)$eigengenes
  else
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=F)$eigengenes

  names(eg.vectors) <- substring(names(eg.vectors),3)
  eg.vectors <- eg.vectors[,centroid.labels]
  return(eg.vectors)
}
getExchangedGenes <- function(old.partition,new.partition){
  stopifnot(length(old.partition) == length(new.partition))
  return(old.partition[old.partition != new.partition])
}

while(exchanged.genes > meg & iteration <= n.iterations){
    print(paste0("Starting partition ",iteration))
    print(paste0("Number of centroids before getting new partition ",ncol(centroids)))
    new.partition <- apply(expr.data,MARGIN=2,getBestModuleCor,centroids=centroids,signed=(net.type == "signed"))
    partitions[[iteration + 1]] <- new.partition
    exchanged.gene.count <- length(getExchangedGenes(partitions[[iteration]], partitions[[iteration + 1]]))
    cat("A total of ", exchanged.gene.count, " genes moved to another partition\n")
    new.partition.in.colors <- centroid.labels[unlist(new.partition)]
    centroids <- getNewCentroids(expr.data,new.partition.in.colors,centroid.labels,min.genes.for.grey)
    exchanged.genes = exchanged.gene.count
    iteration = iteration + 1
}

cat("We finish with",iteration,"iterations\n")


moduleLabels <- partitions[[iteration]]
names(new.partition.in.colors) <- names(moduleLabels)
moduleColors <- new.partition.in.colors
MEs <- centroids

setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clustering/nobatch")
save(moduleLabels, moduleColors, MEs, file = "adjusted_module_info.Rdata")
save(partitions, file = "partitions.Rdata")
################################
#print the modules out to csv  #
################################
dir.create("modules")
setwd("./modules")
colors <- unique(moduleColors)
for (i in 1:length(colors)){
  MyData <- names(which(moduleColors == colors[i]))
  write.csv(MyData, file = paste(colors[i],"_hier.csv", sep = ""))
}
