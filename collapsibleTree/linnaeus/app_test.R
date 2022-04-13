library(shiny)
library(collapsibleTree)
library(data.tree)

#source('linnaeus/app_pancreas.R')
#source('linnaeus/app_neural_larva.R')

if(F){
	master = do_toy_example()
	ct_colors = c("#0000dd","#cd00cd")
	ctypes = c("type1", "type2")
	orit = Clone(master)
	ttt = Clone(master)
	namess = ttt$Get(function(x) if(x$isScar) x$name)
	namess = namess[!is.na(namess)]
	names(namess) = namess
}


do_barplot = function(foc_pie, ct_colors, zeros = FALSE, horiz = FALSE){
 #names(ct_colors) = ctypes
 if(!zeros){
  f_n = foc_pie > 0 
 }else{
  f_n = 1:length(foc_pie)
 }
 foc_pie = foc_pie[f_n]
 #f_n = f_n[order(foc_pie)]
 #foc_pie = foc_pie[f_n]
 large_margin = c(15, 3)
 if(horiz){
	 par(mar=c(rev(large_margin), 0.5, 0.1))
	 barplot(foc_pie, col = ct_colors[f_n], las = 1, space = 0.2, horiz = horiz)
 }else{
	 par(mar=c(large_margin, 0.5, 0.1))
	 barplot(foc_pie, col = ct_colors[f_n], las = 2, space = 0.2)
 }

 grid()
}

# Define server logic required to draw a collapsible tree diagram
cached_trees = list()

get_clicked_node_name = function(dt, input){
	clicked_node = input$node[length(input$node)]
				if(length(input$node) == 0){
					return(dt$name)
				}else{
		   			foc_node = FindNode(dt, clicked_node)
					return(foc_node$name)
				}
}


# Define UI for application that draws a collapsible tree
ui <- fluidPage(
   # Application title
   titlePanel("Linnaeus"),

   # Sidebar with a select input for the root node
   sidebarLayout(position ='right',
      sidebarPanel(
         selectInput("root", "Select a node to render a sub-tree", namess),
         tags$p("Scar from the most recently clicked node:"),
         verbatimTextOutput("str"),
	
	 checkboxInput(inputId = "pieSummary",
	      label = strong("Hide single cells"),  value = TRUE),

	 sliderInput(inputId = "linkLength",
        	label = "Link Length",
	        min = 50, max = 200, value = 150, step = 50),

	 sliderInput(inputId = "treeheight",
        	label = "Tree panel height",
	        min = 500, max = 2000, value = 1000, step = 5)
	),

      # Show a tree diagram with the selected root node
      mainPanel(
        collapsibleTreeOutput("cached", height = '650px'),
	uiOutput("barplot_ctypes.ui")
  	#plotOutput(outputId = "barplot_ctypes", height = "500px")
  	#plotOutput(outputId = "plot.ui", height = "500px")
      )
   )
)

server <- function(input, output){

	observe({
		tree_id <<- paste(input$root, input$pieSummary, input$linkLength)
		if(!(tree_id %in% names(cached_trees))){
			if(input$root != orit$name){
				cached_trees[[tree_id]] <<- collapsibleTree(FindNode(ttt, input$root), collapsed=F, inputId = "node", pieNode=T, pieSummary=input$pieSummary, hide_scars=T, linkLength = input$linkLength, do_collapse= FALSE, ctypes=ctypes, ct_colors=ct_colors)
			} else{
		     	#collapsibleTree(ttt, collapsed=F, inputId = "node", pieNode=T, pieSummary=pieSummary, hide_scars=F, linkLength = linkLength)
			# TODO hide_scars must be set to true as current trees lack the scar field, if F then barplots break as there is no name for a node
		     		cached_trees[[tree_id]] <<- collapsibleTree(ttt, collapsed=F, inputId = "node", pieNode=T, pieSummary=input$pieSummary, hide_scars=T, linkLength = input$linkLength, do_collapse= FALSE, ctypes=ctypes, ct_colors=ct_colors)}
		}
	print(names(cached_trees))
	
	})
	observe(print(names(tree_id)))
	output$cached <- renderCollapsibleTree(cached_trees[[paste(input$root, input$pieSummary, input$linkLength)]])
#	output$plot <- renderCollapsibleTree({			   })

	# Barplot section (and R objects in general)
	observe(
		if(!is.null(input$node)){

			clicked_node = input$node[length(input$node)]
			if(length(input$node) == 0){
   				foc_pie = ttt$pieNode
   				#output$str <- renderPrint(cat(ttt$name))
   				#output$str <- str(as.list(do_summary(ttt)))
	   			foc_node <<- ttt
			}else{
	   			foc_node <<- FindNode(ttt, clicked_node)
	   			foc_pie = foc_node$pieNode
   				#output$str <- renderPrint(cat(foc_node$name))
   				#output$str <- str(as.list(do_summary(foc_node)))
			}
			xSummary = do_summary(foc_node)
			if(sum(foc_pie)>0){
				barplot_width = max(400, sum(foc_pie>0)* 18)
	   			output$barplot_ctypes <- renderPlot(do_barplot(foc_pie, ct_colors = ct_colors))
	   			output$barplot_ctypes.ui <- renderUI({plotOutput(outputId = "barplot_ctypes", height = 600, width = barplot_width)})
			}
   			#output$str <- str(as.list(do_summary(foc_node)))
			output$treeheight <- renderText({input$treeheight})
		}
	)
}

# Run the application
shinyApp(ui = ui, server = server)
