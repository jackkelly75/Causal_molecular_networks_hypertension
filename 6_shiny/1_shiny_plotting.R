library(shiny)
library(igraph)
library(visNetwork)
library(shinythemes)
library(stringr)
library(ppcor)
library(Rgraphviz)
library(RColorBrewer)
library(rsconnect)
library(shiny)

#in dir of the ui.R and server.R
rsconnect::setAccountInfo(name='jack-kelly-manchester', token='7E61F847B764EF2FD9749A736E8FF343', secret='3lyhixJbBcqyCJUeXO7AxUm/IE3Eapi7cXUikF7A')
#runApp()
deployApp(appName = "Causal_molecular_networks_Hypertension")


#run in R
setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/custom_more")
col <- list.dirs(recursive = F)
modules <- vector()
for(i in col){
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/custom_more")
	setwd(i)
	folder <- sub('./', '', i)
	print(paste0("Start ", folder))
	adj <- read.csv("MRPC_adjacency.csv")
	modules <- c(nrow(adj), modules)
}


setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/custom_more")
col <- list.dirs(recursive = F)
modules <- data.frame()
for(i in col){
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/custom_more")
	setwd(i)
	folder <- sub('./', '', i)
	print(paste0("Start ", folder))
	if (file.exists(paste0("MRPC_ancestor_adjacency_", folder, ".csv"))){
		adj <- read.csv(paste0("MRPC_ancestor_adjacency_", folder, ".csv"))
		if(dim(adj)[1] > 1){
			modules[i,1] <- nrow(adj)
			modules[i,2] <- i
		}
	}
}

modules <- modules[modules$V1 > 4,]


#change edge colour by the weight
genkolor <- function(w) {
	k <- vector()
	for (i in 1:length(w)) {
		edge_col <- (1 - w[i])-0.3  #add constant of -0.3 to make the edges a bit darker
		if(edge_col > 1){ #if any of the edges are over 1 then just assign to 1 as needs to be between 0-1
			edge_col <- 1
		}
		k[i] <- rgb(edge_col,edge_col,edge_col)
	}
	return(k)
}

all_nodes = list()
all_edges = list()
num = 1

for(i in modules$V2){
	setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/custom_more")
	setwd(i)
	load("MRPC.fit.Rdata")
	load("MRPC.Rdata")
	folder <- sub('./', '', i)
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
	network_data <- MRPC_data[,ancestors]
	network_partial <- pcor(network_data)$estimate
	rownames(network_partial) <- colnames(network_partial) <- colnames(network_data)
	network_partial[parent_adjacency == 0] <- 0 #get the correlations of the directions
	t.graph <- graph_from_adjacency_matrix(abs(network_partial), weighted = TRUE)
	E(t.graph)$width <- E(t.graph)$weight * 10 #add to weight to scale up to better sizing
	E(t.graph)$color <- "black" #set edge colour to black
	lll <- layout.sugiyama(t.graph)$layout #layout of network
	colors <- brewer.pal(n = 5, name = "Dark2") #get colour pallete to use
	#V(t.graph)$color <- ifelse(V(t.graph)$name == "HYPERTENSION", adjustcolor(colors[2], alpha.f = .8), adjustcolor(colors[1], alpha.f = .8)) #sets color to the node type
	V(t.graph)$color <- ifelse(V(t.graph)$name == "HYPERTENSION", adjustcolor("#bc5090", alpha.f = .8), adjustcolor("#003f5c", alpha.f = .8)) #sets color to the node type
	V(t.graph)$frame.color <- ifelse(V(t.graph)$name == "HYPERTENSION", "#bc5090", "#003f5c")
	#V(t.graph)$frame.color <- ifelse(V(t.graph)$name == "HYPERTENSION", colors[2], colors[1])	
	V(t.graph)$frame.width <- 30 #set size of the border of nodes

	V(t.graph)$size <- 15 #set size of nodes
	#V(t.graph)$color <- ifelse(V(t.graph)$name == "HYPERTENSION", colors[2], colors[1]) #sets color to the node type
	
	#break any that are ENSG into two strings to fit in nodes
	str <- V(t.graph)$name
	str <- gsub("ENSG000", "ENSG000\n", str)
	V(t.graph)$name <- str
	#V(t.graph)$label.cex = 1.1
	
	data <- toVisNetworkData(t.graph)
	nodes <- data$nodes
	edges <- data$edges
	edges$label <- E(t.graph)$width
	#edges$color <- genkolor((E(t.graph)$width)/10)
	edges$color <- "#ffa600"
	edges$weight <- E(t.graph)$weight
	nodes$font.size <- 16 #sets the size of labels in plot
	#nodes$shape= 'circle' #puts text into the node when plotting

	all_nodes[[num]] <- nodes
	all_edges[[num]] <- edges
	names(all_nodes)[num] <- folder
	names(all_edges)[num] <- folder

	num <- num + 1
}

setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/plotting")
saveRDS(all_edges, file = "edges.rds")
saveRDS(all_nodes, file = "nodes.rds")
saveRDS(modules, file = "modules.rds")





setwd("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/MRPC/plotting")
server <- function(input, output) {
	edges <- readRDS(file = "edges.rds")#point this towards data (Rdata file etc.)
	nodes <- readRDS(file = "nodes.rds")#point this towards data (Rdata file etc.)
	darkturquoise <- reactive({
		visNetwork(nodes = nodes$darkturquoise, edges = edges$darkturquoise, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
			visEdges(arrows = "to") %>% 
			#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
 	})
	palevioletred1 <- reactive({
		visNetwork(nodes = nodes$palevioletred1, edges = edges$palevioletred1, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	aliceblue <- reactive({
		visNetwork(nodes = nodes$aliceblue, edges = edges$aliceblue, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	blue <- reactive({
		visNetwork(nodes = nodes$blue, edges = edges$blue, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	coral2 <- reactive({
		visNetwork(nodes = nodes$coral2, edges = edges$coral2, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	darkseagreen2 <- reactive({
		visNetwork(nodes = nodes$darkseagreen2, edges = edges$darkseagreen2, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	firebrick4 <- reactive({
		visNetwork(nodes = nodes$firebrick4, edges = edges$firebrick4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	green <- reactive({
		visNetwork(nodes = nodes$green, edges = edges$green, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	green3 <- reactive({
		visNetwork(nodes = nodes$green3, edges = edges$green3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	indianred4 <- reactive({
		visNetwork(nodes = nodes$indianred4, edges = edges$indianred4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	lightpink4 <- reactive({
		visNetwork(nodes = nodes$lightpink4, edges = edges$lightpink4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	lightskyblue3 <- reactive({
		visNetwork(nodes = nodes$lightskyblue3, edges = edges$lightskyblue3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	pink4 <- reactive({
		visNetwork(nodes = nodes$pink4, edges = edges$pink4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	plum3 <- reactive({
		visNetwork(nodes = nodes$plum3, edges = edges$plum3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	salmon4 <- reactive({
		visNetwork(nodes = nodes$salmon4, edges = edges$salmon4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	skyblue2 <- reactive({
		visNetwork(nodes = nodes$skyblue2, edges = edges$skyblue2, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	tan3 <- reactive({
		visNetwork(nodes = nodes$tan3, edges = edges$tan3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
	})
	graphInput <- reactive({
		switch(input$graph,
		"Aliceblue" = aliceblue(),
		"Blue" = blue(),
		"Coral2" = coral2(),
		"Darkseagreen2" = darkseagreen2(),
		"Darkturquoise" = darkturquoise(),
		"Firebrick4" = firebrick4(),
		"Green" = green(),
		"Green3" = green3(),
		"Indianred4" = indianred4(),
		"Lightpink4" = lightpink4(),
		"Lightskyblue3" = lightskyblue3(),
		"Paleturquoise" = paleturquoise(),
		"Palevioletred1" = palevioletred1(),
		"Pink4" = pink4(),
		"Plum3" = plum3(),
		"Salmon4" = salmon4(),
		"Skyblue2" = skyblue2(),
		"Tan3" = tan3()
		) 
 	})
 	output$selected_graph <- renderVisNetwork({ 
		graphInput()
	})
}

css <- "
#large .selectize-input { line-height: 40px; }
#large .selectize-dropdown { line-height: 30px; }"


ui <- fluidPage(theme = shinytheme("lumen"),
	titlePanel("Causal Networks"),
	#selectInput("graph", "Choose a graph to view:",
    #              choices = c("Darkturquoise", "Palevioletred1")),
    tags$style(type='text/css', css),
    sidebarPanel(width = 2,
	div(id = "large",
    	selectInput("graph", "Choose a graph to view:", list(">5 and <=15 Casusal predictors of Hypertension" = c("Aliceblue", "Blue", "Coral2", "Green3", "Indianred4", "Lightpink4", "Lightskyblue3", "Palevioletred1", "Pink4", "Tan3"), ">15 Casusal predictors of Hypertension" = c("Darkseagreen2", "Darkturquoise", "Firebrick4", "Green", "Paleturquoise", "Plum3", "Salmon4")))
	)
	),

	#visNetworkOutput("selected_graph", height = "800px")
	mainPanel(width = 10,
		visNetworkOutput("selected_graph", height = "800px")
    )
)


shinyApp(ui, server)