library(GeneNet)
library(pcalg) 
library(MRPC)
library(stringr)
library(igraph)
set.seed(120)
library(ppcor)
library(Rgraphviz)
library(psych)


setwd("F:/Manchester-ResearchAssociate/Causal_networks_hypertension/Jan_analysis/eQTL/results/sig")
source("F:/Manchester-ResearchAssociate/Causal_networks_hypertension/Jan_analysis/eQTL/5_MRPC/MRdualPC.R")

col <- list.dirs(recursive = F)
col = col[3:65]

for(module in col){
	setwd("F:/Manchester-ResearchAssociate/Causal_networks_hypertension/Jan_analysis/eQTL/results/sig")
	setwd(module)
	print(module)
	folder <- sub('./', '', module)
	print(folder)
	print("doing MRPC")
	load("MRPC.Rdata")
	p <- length(labels)
	print(paste0("There are ", p, " nodes"))
	GV = sum(str_count(colnames(MRPC_data), pattern = "_b37"))
	skel_pcalg <- ModiSkeleton(data = MRPC_data, suffStat = suffStat_C1, FDR= 0.05, alpha= 0.05, indepTest = "gaussCItest", FDRcontrol = "none", labels = labels, p = p, method = "stable")
	MRPC.fit <- EdgeOrientation(skel_pcalg, GV, suffStat_C1, FDR = 0.05, alpha = 0.05, indepTest = 'gaussCItest', tau = 0.5, lambda = 0.25, FDRcontrol = "none", verbose = FALSE)
	##plotting
	graph <- MRPC.fit@graph
	adjacency <- as(graph, "matrix")
	defAttrs <- getDefaultAttrs()
	remove <- nodes(graph)[-c(GV+1:length(nodes(graph)))]
	temp <- removeNode(remove, graph)
	defAttrs <- getDefaultAttrs()
	igraph <- graph_from_graphnel(temp, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
	d_search = bfs(igraph,"HYPERTENSION",neimode="in", unreachable=FALSE, order=TRUE, dist=TRUE)
	ancestors<-names(d_search$dist[!is.na(d_search$dist)])
	nodes <- nodes(MRPC.fit@graph)
	remove_all <- nodes[!(nodes %in% ancestors)]
	temp <- removeNode(remove_all, MRPC.fit@graph)
	defAttrs <- getDefaultAttrs()
	setwd("F:/Manchester-ResearchAssociate/Causal_networks_hypertension/Jan_analysis/eQTL/results/sig/0_compare_them")
	pdf(paste0("MRPC_", folder, "ancestors.pdf"), height = 30, width = 30) 
	plot (temp, attrs=list(node=list(fontsize=40)))
	dev.off()
	##now do MRdualPC
	print("doing MR dual PC")
	rm(MRPC.fit)
	rm(skel_pcalg)
	filtered <- find_filtered(MRPC_data, 0.05, cor_mat = suffStat_C1$C)
	skel_pcalg <- ModiSkeleton1(data = MRPC_data, filtered = filtered, suffStat = suffStat_C1, FDR= 0.05, alpha= 0.05, indepTest = "gaussCItest", FDRcontrol = "none", labels = labels, p = p, method = "stable")
	MRPC.fit <- EdgeOrientation(skel_pcalg, GV, suffStat_C1, FDR = 0.05, alpha = 0.05, indepTest = 'gaussCItest', tau = 0.5, lambda = 0.25, FDRcontrol = "none", verbose = FALSE)
	##plotting
	graph <- MRPC.fit@graph
	adjacency <- as(graph, "matrix")
	defAttrs <- getDefaultAttrs()
	remove <- nodes(graph)[-c(GV+1:length(nodes(graph)))]
	temp <- removeNode(remove, graph)
	defAttrs <- getDefaultAttrs()
	igraph <- graph_from_graphnel(temp, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
	d_search = bfs(igraph,"HYPERTENSION",neimode="in", unreachable=FALSE, order=TRUE, dist=TRUE)
	ancestors<-names(d_search$dist[!is.na(d_search$dist)])
	nodes <- nodes(MRPC.fit@graph)
	remove_all <- nodes[!(nodes %in% ancestors)]
	temp <- removeNode(remove_all, MRPC.fit@graph)
	defAttrs <- getDefaultAttrs()
	pdf(paste0("MRdualPC_", folder, "ancestors.pdf"), height = 30, width = 30) 
	plot (temp, attrs=list(node=list(fontsize=40)))
	dev.off()
}