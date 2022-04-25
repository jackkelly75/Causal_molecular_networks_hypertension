library(shiny)
library(igraph)
library(visNetwork)
library(shinythemes)

css <- "
#large .selectize-input { line-height: 40px; }
#large .selectize-dropdown { line-height: 30px; }"


ui <- fluidPage(theme = shinytheme("lumen"),
	titlePanel("Causal Networks"),
    tags$style(type='text/css', css),
    sidebarPanel(width = 2,
	div(id = "large",
    	selectInput("graph", "Choose a network to view:", list(">5 and <=15 Casusal predictors of Hypertension" = c("Aliceblue", "Blue", "Coral2", "Green3", "Indianred4", "Lightpink4", "Lightskyblue3", "Palevioletred1", "Pink4", "Tan3"), ">15 Casusal predictors of Hypertension" = c("Darkseagreen2", "Darkturquoise", "Firebrick4", "Green", "Paleturquoise", "Plum3", "Salmon4")))
	),
	# Button
	br(),
	br(),
	p("Download weighted adjacency matrix with partial correlations as .csv files"),
	br(),
	downloadButton("darkturquoise", "Darkturquoise"),
	br(),
	downloadButton("palevioletred1", "Palevioletred1"),
	br(),
	downloadButton("aliceblue", "Aliceblue"),
	br(),
	downloadButton("blue", "Blue"),
	br(),
	downloadButton("coral2", "Coral2"),
	br(),
	downloadButton("darkseagreen2", "Darkseagreen2"),
	br(),
	downloadButton("firebrick4", "Firebrick4"),
	br(),
	downloadButton("green", "Green"),
	br(),
	downloadButton("green3", "Green3"),
	br(),
	downloadButton("indianred4", "Indianred4"),
	br(),
	downloadButton("lightpink4", "Lightpink4"),
	br(),
	downloadButton("lightskyblue3", "Lightskyblue3"),
	br(),
	downloadButton("paleturquoise", "Paleturquoise"),
	br(),
	downloadButton("pink4", "Pink4"),
	br(),
	downloadButton("plum3", "Plum3"),
	br(),
	downloadButton("salmon4", "Salmon4"),
	br(),
	downloadButton("skyblue2", "Skyblue2"),
	br(),
	downloadButton("tan3", "Tan3")
	),
	mainPanel(width = 10,
		visNetworkOutput("selected_graph", height = "1000px")
    )
)

