library(shiny)
library(igraph)
library(visNetwork)
library(shinythemes)

function(input, output) {
	edges <- readRDS(file = "edges.rds")#point this towards data (Rdata file etc.)
	nodes <- readRDS(file = "nodes.rds")#point this towards data (Rdata file etc.)
	all_matrix <- readRDS(file = "all_matrix.rds")
	darkturquoise <- reactive({
		visNetwork(nodes = nodes$darkturquoise, edges = edges$darkturquoise, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
			visEdges(arrows = "to") %>% 
			#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
 	})
	paleturquoise <- reactive({
		visNetwork(nodes = nodes$paleturquoise, edges = edges$paleturquoise, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	palevioletred1 <- reactive({
		visNetwork(nodes = nodes$palevioletred1, edges = edges$palevioletred1, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	aliceblue <- reactive({
		visNetwork(nodes = nodes$aliceblue, edges = edges$aliceblue, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
    		#visOptions(highlightNearest = list(enabled = T, hover = T), manipulation = TRUE)
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	blue <- reactive({
		visNetwork(nodes = nodes$blue, edges = edges$blue, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	coral2 <- reactive({
		visNetwork(nodes = nodes$coral2, edges = edges$coral2, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	darkseagreen2 <- reactive({
		visNetwork(nodes = nodes$darkseagreen2, edges = edges$darkseagreen2, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	firebrick4 <- reactive({
		visNetwork(nodes = nodes$firebrick4, edges = edges$firebrick4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	green <- reactive({
		visNetwork(nodes = nodes$green, edges = edges$green, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	green3 <- reactive({
		visNetwork(nodes = nodes$green3, edges = edges$green3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	indianred4 <- reactive({
		visNetwork(nodes = nodes$indianred4, edges = edges$indianred4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	lightpink4 <- reactive({
		visNetwork(nodes = nodes$lightpink4, edges = edges$lightpink4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	lightskyblue3 <- reactive({
		visNetwork(nodes = nodes$lightskyblue3, edges = edges$lightskyblue3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	pink4 <- reactive({
		visNetwork(nodes = nodes$pink4, edges = edges$pink4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	plum3 <- reactive({
		visNetwork(nodes = nodes$plum3, edges = edges$plum3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	salmon4 <- reactive({
		visNetwork(nodes = nodes$salmon4, edges = edges$salmon4, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	skyblue2 <- reactive({
		visNetwork(nodes = nodes$skyblue2, edges = edges$skyblue2, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
	})
	tan3 <- reactive({
		visNetwork(nodes = nodes$tan3, edges = edges$tan3, height = "4000", width = "100%") %>%
			visIgraphLayout("layout_with_sugiyama") %>%
    		visEdges(arrows = "to", color = "#ffa600") %>% 
    		#visNodes(color = V(t.graph)$color)  %>%
			visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical", degree = list(from = 0, to = 50)), collapse = TRUE, manipulation = TRUE)
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

	output$aliceblue <- downloadHandler(
		filename = function(){paste("aliceblue.csv", sep = "")},
		content = function(file){write.csv(all_matrix$aliceblue, file, row.names = T)}
	)
	output$blue <- downloadHandler(
		filename = function(){paste("blue.csv", sep = "")},
		content = function(file){write.csv(all_matrix$blue, file, row.names = T)}
	)
	output$coral2 <- downloadHandler(
		filename = function(){paste("coral2.csv", sep = "")},
		content = function(file){write.csv(all_matrix$coral2, file, row.names = T)}
	)
	output$darkseagreen2 <- downloadHandler(
		filename = function(){paste("darkseagreen2.csv", sep = "")},
		content = function(file){write.csv(all_matrix$darkseagreen2, file, row.names = T)}
	)
	output$darkturquoise <- downloadHandler(
		filename = function(){paste("darkturquoise.csv", sep = "")},
		content = function(file){write.csv(all_matrix$darkturquoise, file, row.names = T)}
	)
	output$firebrick4 <- downloadHandler(
		filename = function(){paste("firebrick4.csv", sep = "")},
		content = function(file){write.csv(all_matrix$firebrick4, file, row.names = T)}
	)
	output$green <- downloadHandler(
		filename = function(){paste("green.csv", sep = "")},
		content = function(file){write.csv(all_matrix$green, file, row.names = T)}
	)
	output$green3 <- downloadHandler(
		filename = function(){paste("green3.csv", sep = "")},
		content = function(file){write.csv(all_matrix$green3, file, row.names = T)}
	)
	output$indianred4 <- downloadHandler(
		filename = function(){paste("indianred4.csv", sep = "")},
		content = function(file){write.csv(all_matrix$indianred4, file, row.names = T)}
	)
	output$lightpink4 <- downloadHandler(
		filename = function(){paste("lightpink4.csv", sep = "")},
		content = function(file){write.csv(all_matrix$lightpink4, file, row.names = T)}
	)
	output$lightskyblue3 <- downloadHandler(
		filename = function(){paste("lightskyblue3.csv", sep = "")},
		content = function(file){write.csv(all_matrix$lightskyblue3, file, row.names = T)}
	)
	output$paleturquoise <- downloadHandler(
		filename = function(){paste("paleturquoise.csv", sep = "")},
		content = function(file){write.csv(all_matrix$paleturquoise, file, row.names = T)}
	)
	output$pink4 <- downloadHandler(
		filename = function(){paste("pink4.csv", sep = "")},
		content = function(file){write.csv(all_matrix$pink4, file, row.names = T)}
	)
	output$plum3 <- downloadHandler(
		filename = function(){paste("plum3.csv", sep = "")},
		content = function(file){write.csv(all_matrix$plum3, file, row.names = T)}
	)
	output$salmon4 <- downloadHandler(
		filename = function(){paste("salmon4.csv", sep = "")},
		content = function(file){write.csv(all_matrix$salmon4, file, row.names = T)}
	)
	output$skyblue2 <- downloadHandler(
		filename = function(){paste("skyblue2.csv", sep = "")},
		content = function(file){write.csv(all_matrix$skyblue2, file, row.names = T)}
	)
	output$tan3 <- downloadHandler(
		filename = function(){paste("tan3.csv", sep = "")},
		content = function(file){write.csv(all_matrix$tan3, file, row.names = T)}
	)
	output$palevioletred1 <- downloadHandler(
		filename = function(){paste("palevioletred1.csv", sep = "")},
		content = function(file){write.csv(all_matrix$palevioletred1, file, row.names = T)}
	)

}