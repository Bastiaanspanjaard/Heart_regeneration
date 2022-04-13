#' @rdname collapsibleTree
#' @method collapsibleTree Node
#' @export
collapsibleTree.Node <- function(df, hierarchy_attribute = "level",
                                 root = df$name, inputId = NULL, attribute = "leafCount",
                                 aggFun = sum, fill = "lightsteelblue",
                                 linkLength = NULL, fontSize = 10, tooltip = FALSE,
                                 tooltipHtml = NULL,nodeSize = NULL, collapsed = TRUE,
                                 zoomable = TRUE, width = NULL, height = NULL,
				# linnaeus specific
    				nodeSize_sc = 2, nodeLabel_sc = FALSE,
				ct_colors = NULL, ctypes = NULL, sort_by_ctype = TRUE, 
				nodeSize_class = c(   10, 15, 20, 35),
				nodeSize_breaks = c( 0, 5, 20, 500, 1e6),
				angle = 0,
				hide_scars = FALSE,
				pieSummary = TRUE,
				pieNode = FALSE,  
			# linnaeus test
		 	use_scar_as_name = TRUE,
			do_collapse = TRUE,
				...) {

  # acceptable inherent node attributes
  nodeAttr <- c("leafCount", "count")

  # reject bad inputs
  if(!is(df) %in% "Node") stop("df must be a data tree object")
  if(!is.character(fill)) stop("fill must be a either a color or column name")
  if(!is.null(tooltipHtml)) if(!(tooltipHtml %in% df$fields)) stop("tooltipHtml column name is incorrect")
  if(!is.null(nodeSize)) if(!(nodeSize %in% c(df$fields, nodeAttr))) stop("nodeSize column name is incorrect")

  # calculate the right and left margins in pixels
  leftMargin <- nchar(root)
  # This old line is terribly slow
  #rightLabelVector <- df$Get("name", filterFun = function(x) x$level==df$height)
  # Changed to new version 
  rightLabelVector <- Get(df$children, "name", filterFun = function(x) x$level==df$height)
  rightMargin <- max(sapply(rightLabelVector, nchar))

  # Deriving hierarchy variable from data.tree input
  hierarchy <- unique(ToDataFrameTree(df, hierarchy_attribute)[[hierarchy_attribute]])
  if(length(hierarchy) <= 1) stop("hierarchy vector must be greater than length 1")

  # Dealing with ctype colour attribution
  if(is.null(ct_colors)){ct_colors = linnaeus.colors_larva$color}
  if(is.null(ctypes)){ctypes = linnaeus.colors_larva$Cell.type}

# PO this can maybe go away as we are currently using factor() 
  ctypes = ctypes[!is.na(ctypes)]
  ctypes = ctypes[is.na(match(ctypes, "NA"))]

  # create a list that contains the options
  options <- list(
    hierarchy = hierarchy,
    input = inputId,
    attribute = attribute,
    linkLength = linkLength,
    fontSize = fontSize,
    tooltip = tooltip,
    collapsed = collapsed,
    zoomable = zoomable,
    margin = list(
      top = 20,
      bottom = 20,
      left = (leftMargin * fontSize/2) + 25,
      right = (rightMargin * fontSize/2) + 25),
    # linnaeus 
    pieNode = pieNode, 
    useColors = !is.null(ct_colors), 
    colors = ct_colors, 
    angle  = angle,
    do_collapse = do_collapse, # toggle collapsible capabilities
    nodeLabel_sc = ifelse(is.null(nodeLabel_sc), TRUE, nodeLabel_sc )
  )

  # these are the fields that will ultimately end up in the json
  jsonFields <- c("scar")

  if(fill %in% df$fields) {
    # fill in node colors based on column name
    df$Do(function(x) x$fill <- x[[fill]])
    jsonFields <- c(jsonFields, "fill")
  } else {
    # default to using fill value as literal color name
    options$fill <- fill
  }

  # PO determine size classes
  #  disable to keep original tree order
  df = Clone(df)
  if(sort_by_ctype){ 
    SortNumeric(df, decreasing = TRUE, recursive = TRUE, attribute = function(x){
        ifelse(is.na(x$Cell.type) |  x$Cell.type == "NA", "1e4", as.numeric(match(x$Cell.type, ctypes)))
    })
   }

  if(pieNode){
    jsonFields = get_pieNode(df, ctypes = ctypes, nodeSize_breaks = nodeSize_breaks, 
                            nodeSize_sc = nodeSize_sc, jsonFields = jsonFields)
  }



  # linnaeus Test feature HideScarnames -> show them on tooltip
  if(hide_scars){    
    df = Clone(df)
    tra  = data.tree::Traverse(df, 'level')
    sapply(tra, function(x){
# TODO determine a fixed expected node field list e.g. scar, name, is, Scar, etc.
      x$scar = gsub("_", ".", x$name)
      xSummary = do_summary(x)
      tot_scells = ifelse(!is.null(x$pieNode), sum(x$pieNode), xSummary$progeny)  
      x$tp <- sprintf('<h4 style="color: #2e6c80;">%s</h4><br><h4 style="color: #2e6c80;">shown:\t\t%s<br>in node:\t\t%s</h4>', 
          ifelse(x$isScar, x$scar, x$Cell.type), tot_scells, xSummary$progeny )
      x$scar = x$name}
    )
    options$tooltip = TRUE
    tooltip = TRUE
    tooltipHtml = "tp"
  }

  if(pieSummary & pieNode){
   # Only after collecting the statistics for the scar nodes we get rid of the scells
    t <- data.tree::Traverse(df, 'post-order')
   	    data.tree::Do(t, function(x) {
            if(x$isLeaf & !x$isRoot & !x$isScar){	
                x$parent$RemoveChild(x$name)
            }
    })  
  }
  if("Main" %in% df$fieldsAll){

    jsonFields <- c(jsonFields, "Main")
  
  }

  # only necessary to perform these calculations if there is a tooltip
  if(tooltip & is.null(tooltipHtml)) {
    t <- data.tree::Traverse(df, hierarchy_attribute)
    if(substitute(identity)=="identity") {
      # for identity, leave the tooltips as is
      data.tree::Do(t, function(x) {
        x$WeightOfNode <- x[[attribute]]
      })
    } else {
      # traverse down the tree and compute the weights of each node for the tooltip
      data.tree::Do(t, function(x) {
        x$WeightOfNode <- data.tree::Aggregate(x, attribute, aggFun)
        # make the tooltips look nice
        x$WeightOfNode <- prettyNum(
          x$WeightOfNode, big.mark = ",", digits = 3, scientific = FALSE
        )
      })
    }
    jsonFields <- c(jsonFields, "WeightOfNode")
  }

  # if tooltipHtml is specified, pass it on in the data
  if(tooltip & !is.null(tooltipHtml)) {
    df$Do(function(x) x$tooltip <- x[[tooltipHtml]])
    jsonFields <- c(jsonFields, "tooltip")
  }

  # only necessary to perform these calculations if there is a nodeSize specified
  if(!is.null(nodeSize) & !pieNode) {
    # Scale factor to keep the median leaf size around 10
    scaleFactor <- 10/data.tree::Aggregate(df, nodeSize, stats::median)
    t <- data.tree::Traverse(df, hierarchy_attribute)
    # traverse down the tree and compute the size of each node
    data.tree::Do(t, function(x) {
      # x$SizeOfNode <- data.tree::Aggregate(x, nodeSize, aggFun)
      # scale node growth to area rather than radius and round
      # x$SizeOfNode <- round(sqrt(x$SizeOfNode*scaleFactor)*pi, 2)
      x$SizeOfNode <- eval(parse(text = paste("x$", nodeSize, sep = "")))
      # scale node growth to area rather than radius and round
      x$SizeOfNode <- round(sqrt(x$SizeOfNode)*pi, 2)
    })
    # update left margin based on new root size
    jsonFields <- c(jsonFields, "SizeOfNode")
  }

  # keep only the JSON fields that are necessary
  if(is.null(jsonFields)) jsonFields <- NA
  data <- data.tree::ToListExplicit(df, unname = TRUE, keepOnly = jsonFields)
  if(use_scar_as_name) {data <- renameNode(data)}
  # pass the data and options using 'x'
  x <- list(
    data = data,
    options = options
  )

  # create the widget
  htmlwidgets::createWidget(
    "collapsibleTree", x, width = width, height = height,
    htmlwidgets::sizingPolicy(viewer.padding = 0)
  )
}

