# Linnaeus helper functios

#' Helper function for pieNode modality of a Linnaeus data.tree object
#' @param df data.tree object to compute cell type counts. Requires Cell.type field
#' @param ctypes cell types
#' @param nodeSize_breaks pie node size breaks
#' @param nodeSize_sc node size for single cells (currently overriden by JS)
#' @param jsonFields what to JSON
#' @rdname collapsibleTree
#' @export
get_pieNode = function(df, ctypes, 
				nodeSize_class = c(   10, 15, 20, 35), 
				nodeSize_breaks = c( 0, 5, 20, 100, 1e6), 
				nodeSize_sc = 2, jsonFields  = NULL
			){
    t <- data.tree::Traverse(df, 'level')
    data.tree::Do(t, function(x) {
	x$isScar = !x$isLeaf & !x$isRoot & x$Cell.type == "NA" # This is a source of problems
	if(x$isRoot) {
		x$isScar = TRUE
		x$Cell.type = "_"
	}
	xpieNode = x$Get("Cell.type")
	x$ct = x$Cell.type
	x$pieNode = table(factor(array(xpieNode), levels=ctypes))
	x$SizeOfNode = nodeSize_class[cut(sum(x$pieNode), breaks=nodeSize_breaks, include.lowest=T, labels=F )]
	if(!x$isScar) {
		x$SizeOfNode = nodeSize_sc
	}else{
		sapply(x$children, function(child){child$parSize = length(x$children)}) 
	}
    })

    if(!is.null(jsonFields)){
      jsonFields <- c(jsonFields, "pieNode")
      jsonFields <- c(jsonFields, "parSize") # keeps a record of the size of parent; used to decide wheter to show cell type of single cell
      jsonFields <- c(jsonFields, "ct")
      jsonFields <- c(jsonFields, "SizeOfNode")
      jsonFields <- c(jsonFields, "isScar")
      return(jsonFields)
    }
}



#' @rdname collapsibleTree
#' @export
pieProportions <- function(node) {
  return(c(node$Cell.type, sapply(node$children, pieProportions)))
}

#' @rdname collapsibleTree
#' @export
# color legend
do_color_map <- function(name, ct_colors, ctypes, width = 2.5, height=10){
	pdf(paste0(name,'_color_map.pdf'), width=width, height=height)
	par(mar=c(1.1, 10, 1.1, 1.1))
	image(y=1:length(ct_colors), x=1, t(as.matrix(1:length(ct_colors))), col = ct_colors, axes=F, ylab='', xlab='') 
		axis(2, at=1:length(ct_colors), las=2, labels = ctypes, cex.axis=0.5)
	dev.off()
}
#' Helper function to sort children by attribute, casting the given value to numeric
#' @rdname collapsibleTree
#' @export
SortNumeric = function (node, attribute, ..., decreasing = FALSE, recursive = TRUE)
{
    if (node$isLeaf)
        return()
    ChildL <- sapply(node$children, function(x) GetAttribute(x,
        attribute, ...))
    names(ChildL) <- names(node$children)
    node$children <- node$children[order(as.numeric(ChildL), decreasing = decreasing,
        na.last = TRUE)]
    if (recursive)
        for (child in node$children) SortNumeric(child, attribute, ...,
            decreasing = decreasing, recursive = recursive)
    invisible(node)
}


#' @rdname collapsibleTree
#' @export
renameNode <- function(node){
  if("scar" %in% names(node)){
    node$name <- node$scar
  }else{
    node$name <- ""
  }

  if("children" %in% names(node)){
    for(i in 1:length(node$children)){
      node$children[[i]] <- renameNode(node$children[[i]])
    }
  }

  return(node)
}

#' Node Summary
#' @rdname collapsibleTree
#' @export 
do_summary = function(x){
	if(!is.null(x$pieNode)){
		# the last -1 is to not count for x itself
		progeny = sum(x$Get(function(x){if(!x$isLeaf){length(x$children)-1}}), na.rm=T) +1
		name = x$name
	return(list(progeny = progeny,
		name = name
	))
	}	
}

#' Helper function for pieNode modality of a Linnaeus data.tree object
#' @rdname collapsibleTree
#' @export
linnaeus.sets = function(setname='adult'){
	x = system.file("extdata", linnaeus.sets[[name]], package = "collapsibleTree")
	x = read.table(x, header = T, stringsAsFactors = FALSE, sep = '\t', comment.char='')
	return(x)

#write.table(color_larvea, file='src/linnaeus-scripts/collapsibleTree/inst/extdata/colors_larvae.txt', quote=F, sep='\t', row.names=F)
}

#' Function to generate a list containing cached instances of different trees and parameters.
#' @param df data.tree object to compute cell type counts. Requires Cell.type field
#' @rdname collapsibleTree
#' @export 



#' @rdname collapsibleTree
#' @export 
do_toy_example = function(type1="type1", type2="type2", scar_color="#aaaaaa", type1_color = "#ff0000", type2_color="#0000ff", scarSize=6, scSize = 2.5, width=500, height=500, return_widget = FALSE){
master <- Node$new('master', size=scarSize, fill=scar_color, Cell.type="NA")
 scarA <- master$AddChild('scarA', size=scarSize, fill=scar_color, Cell.type='NA')
	scarC <-  scarA$AddChild('scarC', size=scarSize, fill=scar_color, Cell.type='NA')
			scC1 = scarC$AddChild('scC1', size=scSize, fill=type1_color, Cell.type=type1)
			scC2 = scarC$AddChild('scC2', size=scSize, fill=type1_color, Cell.type=type1)
			scC3 = scarC$AddChild('scC3', size=scSize, fill=type1_color, Cell.type=type1)
			scC4 = scarC$AddChild('scC4', size=scSize, fill=type1_color, Cell.type=type1)
	scA1 = scarA$AddChild('scA1', size=scSize, fill=type1_color, Cell.type=type1)
	scA2 = scarA$AddChild('scA2', size=scSize, fill=type1_color, Cell.type=type1)
	scA3 = scarA$AddChild('scA3', size=scSize, fill=type2_color, Cell.type=type2)
	scA4 = scarA$AddChild('scA4', size=scSize, fill=type2_color, Cell.type=type2)
	scarD <-  scarA$AddChild('scarD', size=scarSize, fill=scar_color, Cell.type='NA')
			scD1 = scarD$AddChild('scD1', size=scSize, fill=type1_color, Cell.type=type1)
			scD2 = scarD$AddChild('scD2', size=scSize, fill=type1_color, Cell.type=type1)
			scD3 = scarD$AddChild('scD3', size=scSize, fill=type2_color, Cell.type=type2)
			scD4 = scarD$AddChild('scD4', size=scSize, fill=type2_color, Cell.type=type2)
 sc1 = master$AddChild('sc1', size=scSize, fill=type1_color, Cell.type=type1)
 sc2 = master$AddChild('sc2', size=scSize, fill=type1_color, Cell.type=type1)
 sc3 = master$AddChild('sc3', size=scSize, fill=type2_color, Cell.type=type2)
 sc4 = master$AddChild('sc4', size=scSize, fill=type2_color, Cell.type=type2)
 scarB <- master$AddChild('scarB', size=scarSize, fill=scar_color, Cell.type='NA')
	scarE <-  scarB$AddChild('scarE', size=scarSize, fill=scar_color, Cell.type='NA')
			scE1 = scarE$AddChild('scE1', size=scSize, fill=type1_color, Cell.type=type1)
			scE2 = scarE$AddChild('scE2', size=scSize, fill=type1_color, Cell.type=type1)
			scE3 = scarE$AddChild('scE3', size=scSize, fill=type2_color, Cell.type=type2)
			scE4 = scarE$AddChild('scE4', size=scSize, fill=type2_color, Cell.type=type2)
	scB1 = scarB$AddChild('scB1', size=scSize, fill=type1_color, Cell.type=type1)
	scB2 = scarB$AddChild('scB2', size=scSize, fill=type1_color, Cell.type=type1)
	scB3 = scarB$AddChild('scB3', size=scSize, fill=type2_color, Cell.type=type2)
	scB4 = scarB$AddChild('scB4', size=scSize, fill=type2_color, Cell.type=type2)
	scarF <-  scarB$AddChild('scarF', size=scarSize, fill=scar_color, Cell.type='NA')
			scF1 = scarF$AddChild('scF1', size=scSize, fill=type2_color, Cell.type=type2)
			scF2 = scarF$AddChild('scF2', size=scSize, fill=type2_color, Cell.type=type2)
			scF3 = scarF$AddChild('scF3', size=scSize, fill=type2_color, Cell.type=type2)
			scF4 = scarF$AddChild('scF4', size=scSize, fill=type2_color, Cell.type=type2)
# TODO sort children to have symmetrical trees
	get_pieNode(master, ctypes=c(type1, type2))
	if(return_widget){
		widget_sc = collapsibleTree(master, collapsed = F, pieSummary=F, pieNode=F, ctypes=unique(master$Get("Cell.type")), nodeSize='size', fill="fill", ct_colors=c(type2_color,type1_color), width=width, height=height, sort_by_ctype = FALSE)
		widget_summary = collapsibleTree(master, collapsed = F, pieSummary=T, pieNode=T, ctypes=unique(master$Get("Cell.type")), ct_colors=c(type1_color,type2_color), width=width, height=height)
		return(list(widget_sc, widget_summary))
	}else{
		return(master)
	}
}
