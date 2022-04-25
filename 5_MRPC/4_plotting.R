library(stringr)
library(Rgraphviz)
library(igraph)
library(RColorBrewer)
library(ppcor)

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

for(i in modules){
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/nobatch")
	setwd(i)
	load("MRPC.fit.Rdata")
	load("MRPC.Rdata")
	defAttrs <- getDefaultAttrs()
	GV = sum(str_count(colnames(MRPC_data), pattern = "_b37"))
	graph <- MRPC.fit@graph
	adjacency <- as(graph, "matrix")
	remove <- nodes(graph)[-c(GV+1:length(nodes(graph)))]
	temp <- removeNode(remove, graph)

	igraph <- graph_from_graphnel(temp, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
	d_search = bfs(igraph,"HYPERTENSION",neimode="in", unreachable=FALSE, order=TRUE, dist=TRUE)
	ancestors<-names(d_search$dist[!is.na(d_search$dist)])
	nodes <- nodes(MRPC.fit@graph)
	remove_all <- nodes[!(nodes %in% ancestors)]
	temp <- removeNode(remove_all, MRPC.fit@graph)
	parent_adjacency <- adjacency[ancestors,ancestors] #gets the directional adjacency matrix wth 0 where no edge (used for plotting)
	#get edge weights by marginal correlation
	#parent_corr <- suffStat_C1$C[ancestors,ancestors]
	#parent_corr[parent_adjacency == 0] <- 0 #get the correlations of the directions
	#write.csv(parent_corr, quote= F , file = "MRPC_weighted_edges.csv")
	#get edge weights by partial correlation
	network_data <- MRPC_data[,ancestors]
	network_partial <- pcor(network_data)$estimate
	rownames(network_partial) <- colnames(network_partial) <- colnames(network_data)
	network_partial[parent_adjacency == 0] <- 0 #get the correlations of the directions

	#t.graph <- graph_from_adjacency_matrix(abs(parent_corr), weighted = TRUE)
	t.graph <- graph_from_adjacency_matrix(abs(network_partial), weighted = TRUE)
	E(t.graph)$weight
	E(t.graph)$width <- E(t.graph)$weight * 10 #add to weight to scale up to better sizing
	#info on layouts - https://igraph.org/r/doc/layout_.html
	lll <- layout.sugiyama(t.graph)$layout
	#lll <- layout_with_sugiyama(t.graph, hgap = 5, vgap = 5)$layout

	#colourBrewer palletes - https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
	colors <- brewer.pal(n = 5, name = "Dark2")
	V(t.graph)$color <- ifelse(V(t.graph)$name == "HYPERTENSION", colors[2], colors[1])
	#V(t.graph)$label.cex = 2

	#change size based on the number of nodes
	#length(V(t.graph))
	#if smaller than cutoff (set 3) then sclae up size and tezxt size

	str <- V(t.graph)$name
	str <- gsub("ENSG000", "ENSG000\n", str)
	V(t.graph)$name <- str

	V(t.graph)$label.cex = 1.1
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/nobatch")
	folder <- sub('./', '', i)
	png(paste0("MRPC_", folder,"_image.png"), width = 10, height =10 , units = "in", res = 800) 
	plot(t.graph, layout = lll, edge.color=colors[3], vertex.frame.color="#ffffff", vertex.label.color="black")
	dev.off()
}