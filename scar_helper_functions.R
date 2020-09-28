# Description ####
# Functions for analysis of scars and WT

# Written by B. Spanjaard, 2017

# Dependencies ####
suppressPackageStartupMessages(library(ggplot2))
# require(proxy)
# require(tsne)
# require(dbscan)
require(igraph)
suppressPackageStartupMessages(library(reshape2))
# require(scales)
require(stringdist)
require(ape)
require(poisbinom)
require(plyr)
require(data.tree)
require(treemap)
require(collapsibleTree)
require(data.table)
# source("./Scripts/linnaeus-scripts/collapsibleTree.R")
# source("./Scripts/linnaeus-scripts/collapsibleTree.data.tree.R")


# Parameters ####
wildtype.seq <- 
  "GTTCAAGACCATCTACATGGCCAAGAAGCCCGTGCAACTGCCCGGCTACTACTACGTGGACACCAAGCTGGACAT"
sc.primer <- "GTTCAAGACCATCTACATGGCC"
sc.primer.length <- nchar(sc.primer)
wildtype.clip <- substring(wildtype.seq, sc.primer.length + 1, nchar(wildtype.seq) - 2)
theme_bs <- 
  theme_update(text = element_text(size = 24),
               panel.background = element_rect(fill = "white", colour = "black"),
               panel.grid.major.y = element_line(colour = "grey"),
               panel.grid.major.x = element_line(colour = "grey"),
               legend.key = element_rect(fill = "white"),
               plot.title = element_text(hjust = 0.5))

# Functions ####
binary <- function(x) {
  i <- 0
  string <- numeric(32)
  while(x > 0) {
    string[32 - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  first <- match(1, string)
  string[first:32] 
}

calculate.adjacency <- function(scars, scar.probabilities, wt.seq = wildtype.seq,
                                wt.remove = T){
  # Calculate the adjacency matrix of lineage connections between all cells
  # with scars in [scars], based on their scars (in [scars]) and the weight of 
  # those scars in [scar.weights].
  
  # The object [scars] has a row per scar in a cell, and at least the following
  # columns: cell names ($Cell), CIGAR codes ($CIGAR), locations ($Location),
  # and sequence ($Sequence). The object [scar.weights] has a row per scar,
  # and at least the following columns: sequence ($Sequence), and p-value 
  # for the creation of that scar ($p).
  
  # Append the probabilities to all scars found, remove wildtype and find the
  # p-value of the connections between the cells, based on the scars they
  # have and the probabilities of those scars. The p-value of a connection
  # is the product of the probabilities of the constituent scars, reflecting
  # the probability that the observed connection is due to scars being
  # independently created in two cells.
  # 
  # 
  # named.scars <- get.scar.names(scars, wt.remove = wt.remove, 
  #                               return.seq.and.loc = T, wt.seq = wt.seq)
  # named.scars <- merge(named.scars, scar.probabilities)
  # connections.p <- find.all.p.links(named.scars)
  # 
  connections.p <- find.all.p.links(scars)
  
  # Add all connections of cells with themselves (something 
  # "find.weighted.links" by design does not return), with p-value 0.
  cells <- unique(scars$Cell)
  self.connections <- data.frame(matrix(nrow = length(cells), ncol = 3))
  colnames(self.connections) <- colnames(connections.p)
  self.connections$Cell.1 <- cells
  self.connections$Cell.2 <- cells
  self.connections$Strength <- 0
  connections.p <- rbind(connections.p, self.connections)
  
  # Create a full adjacency matrix from these connections by casting
  # the resulting list of connections into a NxN-matrix, N the number of cells,
  # containing the p-value of all cell-cell connections.
  adjacency <- dcast(connections.p, Cell.1 ~ Cell.2, value.var = "Strength",
                     fill = 1)
  rownames(adjacency) <- adjacency$Cell.1
  adjacency <- adjacency[, -1]
  
  return(adjacency)
}

calculate.cell.scar.distance <- function(scar.cell.matrix, scar.inner.product.matrix){
  # Calculate the distance between cells based on the cells they have.
  # Distance between two cells:
  # D(c1,c2) = ||c1 - c2|| = sqrt(<c1,c1> - 2<c1,c2> + <c2,c2>), with <c1,c2>
  # similarity/inner product between cells. Let i,j be cell indices, and a,b
  # scar indices, then the similarity matrix between all cells is
  # S_ij = C_iaI_abC_bj, with C_ai the scar-cell matrix, and I_ab the scar
  # inner product matrix.
  # Ingredients for distances that are accepted: 
  # 1) scar-cell matrix 
  #   a) observed (binary values 0 and 1)
  #   b) imputed (values between 0 - absolutely not present and 1 - observed)
  # 2) scar inner product matrix
  #   a) Euclidean
  #   b) Weighted (by scar creation probabilities)
  
  scars.selected <- rownames(scar.cell.matrix)
  inner.product.selected <- 
    scar.inner.product.matrix[scars.selected, scars.selected]
  
  cell.sim <- t(scar.cell.matrix) %*% inner.product.selected %*% scar.cell.matrix
  sim.1 <- matrix(rep(diag(cell.sim), nrow(cell.sim)), nrow = nrow(cell.sim))
  cell.dist.input <- sqrt(sim.1 + t(sim.1) - 2 * cell.sim)
  cell.dist <- as.dist(cell.dist.input)
  
  return(cell.dist)
}

calculate.links <- function(fscars.rcc){
  # To calculate the number of links between cell properties (such as clusters
  # or fish), we first need to establish which cells are linked. Since we know 
  # which cell has which property, we can then easily determine the number of 
  # links between different properties. Establishing which cells are linked is 
  # very much like matrix multiplication of the count matrix of all scars in 
  # all cells with its transpose.
  property <- colnames(fscars.rcc)[3]
  
  cell.scar <- dcast(fscars.rcc, Scar ~ Cell, fun.aggregate = length,
                     value.var = property)
  cell.scar <- cell.scar[complete.cases(cell.scar), ]
  rownames(cell.scar) <- cell.scar$Scar
  cell.scar <- as.matrix(cell.scar[, -1])
  cell.links <- t(cell.scar) %*% cell.scar
  
  # So far, this has been matrix multiplication, but links are a binary entity
  # and links between a cell and itself aren't possible. So we remove the 
  # diagonal and binarize the rest.
  diag(cell.links) <- 0
  cell.links[cell.links > 0] <- 1
  # Plot to check
  # pheatmap::pheatmap(cell.links)
  
  # To avoid overcounting, we remove the upper triangle of the cell links.
  cell.links.lt <- cell.links
  cell.links.lt[upper.tri(cell.links.lt)] <- 0
  
  # And now we convert the cell links to property links by multiplying with the
  # the matrix that designates which cell has which property (and its
  # transpose). We first need to create that matrix though.
  cell.property.long <- unique(fscars.rcc[, c("Cell", property)])
  colnames(cell.property.long)[2] <- "Property"
  cell.property <- dcast(cell.property.long, Cell ~ Property, 
                         fun.aggregate = length, value.var = "Property")
  rownames(cell.property) <- cell.property$Cell
  cell.property <- as.matrix(cell.property[, -1])
  # And once the cell-property matrix is created, we can easily count the number
  # of links between properties. We normalize by the number of links, add the 
  # upper-triangle elements for plotting convenience, and reshape the property
  # links to long format so that we can easily compare these values to the
  # values we get with resampling.
  property.links.lt <- t(cell.property) %*% cell.links.lt %*% cell.property
  property.links.lt.n <- property.links.lt/sum(property.links.lt)
  property.links <- property.links.lt.n
  property.links[upper.tri(property.links)] <- 
    property.links[lower.tri(property.links)]
  property.links.l <- melt(property.links, 
                           varnames = paste(property, 1:2, sep = "."),
                           value.name = "Link.ratio")
  
  return(property.links.l)
}

calculate.links.bos <- function(fscars.rcc){
  # To calculate the number of links between strata with different properties
  # (such as clusters or fish) with the bag of scars technique, we assign every
  # scar to a stratum (i.e. decouple them from the cells they're found in), and
  # calculate the links between and within strata.
  property <- colnames(fscars.rcc)[3]
  
  scars.strata <- fscars.rcc[, -1]
  scars.strata$Scar.container <- 1:nrow(scars.strata)
  scars.container.matrix <- dcast(scars.strata, Scar ~ Scar.container,
                                  fun.aggregate = length,
                                  value.var = property)
  scars.container.matrix <- as.matrix(scars.container.matrix[, -1])
  
  scar.links <- t(scars.container.matrix) %*% scars.container.matrix
  scar.links[upper.tri(scar.links, diag = T)] <- 0
  
  cell.property.long <- unique(scars.strata[, c(property, "Scar.container")])
  colnames(cell.property.long)[1] <- "Property"
  cell.property <- dcast(cell.property.long, Scar.container ~ Property, 
                         fun.aggregate = length, value.var = "Property")
  cell.property <- as.matrix(cell.property[, -1])
  
  # So far, this has been matrix multiplication, but links are a binary entity
  # and links between a cell and itself aren't possible. So we remove the 
  # diagonal and binarize the rest.
  
  #   # To avoid overcounting, we remove the upper triangle of the cell links.
  
  # And now we convert the cell links to property links by multiplying with the
  # the matrix that designates which cell has which property (and its
  # transpose). We first need to create that matrix though.
  
  # And once the cell-property matrix is created, we can easily count the number
  # of links between properties. We normalize by the number of links, add the 
  # upper-triangle elements for plotting convenience, and reshape the property
  # links to long format so that we can easily compare these values to the
  # values we get with resampling.
  property.links.lt <- t(cell.property) %*% scar.links %*% cell.property
  # property.links.lt.n <- property.links.lt/sum(property.links.lt)
  # property.links <- property.links.lt.n
  property.links <- property.links.lt
  property.links[upper.tri(property.links)] <- 
    property.links.lt[upper.tri(property.links.lt)] +
    property.links.lt[lower.tri(property.links.lt)]
  property.links[lower.tri(property.links)] <- 
    property.links[upper.tri(property.links)]
  property.links[upper.tri(property.links)] <- 
    property.links[lower.tri(property.links)]
  property.links.l <- melt(property.links, 
                           varnames = paste(property, 1:2, sep = "."),
                           value.name = "Links")
  
  return(property.links.l)
}

calculate.node.enrichment <- function(node.counts, node.counts.agg,
                                      test.alternative = "two.sided", p.adj.method = "fdr"){
  # Calculate the enrichment of cell types in nodes, using a dataframe of cell-
  # type specific cumulative node counts and a dataframe of aggregated
  # cumulative node counts (the latter one is strictly speaking unnecessary
  # since it could also be calculated from the former).
  # The test done is a binomial test: can the number of cells in a
  # node be explained by random sampling from the parent node? The p-values are
  # corrected for multiple testing.
  # Testing and p-value correcting can be done using all implemented binom.test
  # methods (two-sided, greater, less), and using 
  e.calc <- 
    merge(node.counts[, c("Node", "Cell.type", "Freq")],
          node.counts.agg[, c("Node", "Freq")],
          by = "Node")
  colnames(e.calc)[3:4] <- c("Node.tc", "Node.total")
  e.calc$Parent.node <-
    sapply(e.calc$Node,
           function(x) {
             if(x == "0_1"){
               return(NA)
             }else{
               y <- unlist(strsplit(x, "_"))
               return(paste(y[-length(y)], collapse = "_"))
             }
           }
    )
  e.calc <- e.calc[complete.cases(e.calc), ]
  e.calc <- merge(e.calc, node.counts[, c("Node", "Cell.type", "Freq")],
                  by.x = c("Parent.node", "Cell.type"), by.y = c("Node", "Cell.type"))
  colnames(e.calc)[6] <- "Parent.tc"
  e.calc <- merge(e.calc, node.counts.agg[, c("Node", "Freq")],
                  by.x = "Parent.node", by.y = "Node")
  colnames(e.calc)[7] <- "Parent.total"
  e.calc <- e.calc[, c("Cell.type", "Node", "Parent.node", "Node.tc", "Node.total",
                       "Parent.tc", "Parent.total")]
  e.calc$Binom.p <-
    apply(e.calc[, c("Node.tc", "Node.total", "Parent.tc", "Parent.total")], 1,
          function(x){
            bt <- binom.test(x = x[1], n = x[2], p = x[3]/x[4], 
                             alternative = test.alternative)
            return(bt$p.value)
          }
    )
  e.calc$p.adj <- p.adjust(e.calc$Binom.p, method = p.adj.method)
  
  return(e.calc)
}

connections.for.graph <- function(current.cs){
  # First create a scar-scar matrix that counts how often each combination
  # of scars has been observed.
  scars.in.tree <- data.frame(table(current.cs$Scar))
  colnames(scars.in.tree) <- c("Scar", "Count")
  
  
  # OLD
  # scar.cell <- acast(current.cs, Cell ~ Scar, fun.aggregate = length,
  #                    value.var = "Scar")

  # OLD
  # current.cs$Presence <- 1
  # scar.cell <- dcast(setDF(current.cs), Cell ~ Scar, value.var = "Presence")
  # rownames(scar.cell) <- scar.cell$Cell
  # scar.cell <- scar.cell[, -1]
  # scar.cell[is.na(scar.cell)] <- 0
  # scar.cell <- as.matrix(scar.cell)
  # scar.links <- t(scar.cell) %*% scar.cell
  # Now create a list of all connections between scars.
  # scar.connections <- data.frame(melt(scar.links,
  #                                     varnames = c("Scar.A", "Scar.B"),
  #                                     value.name = "Link"))
  # colnames(scar.connections)[3] <- "x_AB"
  
  # NEW
  # timer.start <- Sys.time()
  current.cs.2 <- current.cs
  colnames(current.cs.2)[2] <- "Scar.B"
  connections.pre <- merge(current.cs[, c("Cell", "Scar")], 
                       current.cs.2[, c("Cell", "Scar.B")], all= T)
  scar.connections <- ddply(connections.pre[, c("Scar", "Scar.B")],
                         .(Scar, Scar.B),nrow)
  colnames(scar.connections) <- c("Scar.A", "Scar.B", "x_AB")
  # timer.end <- Sys.time()
  # timer.time <- timer.end - timer.start
  # print(timer.time)
  
  scar.counts <- 
    scar.connections[scar.connections$Scar.A == scar.connections$Scar.B, 
                     c("Scar.A", "x_AB")]
  colnames(scar.counts) <- c("Scar", "Freq")
  scar.connections <-
    scar.connections[scar.connections$Scar.A != scar.connections$Scar.B, ]
  
  scar.connections <-
    merge(scar.connections, scar.counts, by.x = "Scar.A", by.y = "Scar")
  colnames(scar.connections)[4] <- "x_A"
  scar.connections <-
    merge(scar.connections, scar.counts, by.x = "Scar.B", by.y = "Scar")
  colnames(scar.connections)[5] <- "x_B"
  
  scar.connections <- 
    scar.connections[, c("Scar.A", "Scar.B", "x_A", "x_B", "x_AB")]
  
  return(scar.connections)
}

count.type.connections <- function(connections.type){
  # Count data in long format. This is not symmetric because links only occur
  # once in [connections.type]
  count.table <-
    data.frame(table(connections.type[, c("Cell.type.1", "Cell.type.2")]))
  count.matrix <-
    acast(count.table, Cell.type.1 ~ Cell.type.2, value.var = "Freq")
  symm.count.matrix <- count.matrix
  # To make this matrix symmetric, first add the upper triangle to the lower
  # triangle, then copy the lower triangle to the upper triangle.
  symm.count.matrix[lower.tri(symm.count.matrix)] <-
    symm.count.matrix[lower.tri(symm.count.matrix)] + 
    t(symm.count.matrix)[lower.tri(t(symm.count.matrix))]
  symm.count.matrix[upper.tri(symm.count.matrix)] <- 
    t(symm.count.matrix)[upper.tri(t(symm.count.matrix))]
  
  symm.counts.long <- melt(symm.count.matrix, 
                           varnames = c("Cell.type.1", "Cell.type.2"),
                           value.name = "Count")
  
  return(symm.counts.long)
}

# create.degree.lls(current.cs.component, current.graph)
# cs <- current.cs.component
# graph <- current.graph
create.degree.lls <- function(cs, graph){
  # Create a dataframe that combines count and degree information on all
  # scars, and tests the congruence between the two for each scar, assuming
  # it is the earliest created scar of all scars present.
  
  # Determine scar counts and scar degrees for scars in the current graph.
  graph.degrees <- 
    data.frame(Degree = degree(graph, mode = "all", loops = F))
  graph.degrees$Scar <- rownames(graph.degrees)
  graph.degrees <- graph.degrees[order(-graph.degrees$Degree), ]
  
  # Put together degree and count for all scars
  # scar.counts <- data.frame(table(cs$Scar, cs$Cell.type))
  # colnames(scar.counts) <- c("Scar", "Cell.type", "Scar.count")
  # scar.lls <- graph.degrees
  scar.counts <- data.frame(table(cs$Scar))
  colnames(scar.counts) <- c("Scar", "Scar.count")
  scar.lls <- merge(graph.degrees, scar.counts)

  # Determine scar detection rates, assuming for every scar that it is the first
  # scar created. For p(A) the detection rate of A, we calculate 
  # p(A) = #(cell that have A and another scar)/#(cells that have another scar).
  # This formula is explained in the supplementary methods. The way we calculate
  # these numbers is by using #(A+other) = #(A) - #(only A), and 
  # #(other) = #(A+other) + #(only other). These can be calculated from three
  # elements: #total is the total number of cells with a scar, scar.cell.matrix
  # is the matrix of scars(rows)xcells(columns), and scar.cell1.matrix is similar
  # to scar.cell.matrix but only contains the cells that have only one scar.
  # Then #(A) is rowSums(scar.cell.matrix), #(only A) is rowSums(scar.cell1.matrix),
  # #(only other) is #total - #A.
  scar.cell.matrix <- 
    acast(cbind(cs, 
                data.frame(Presence = rep(1, nrow(cs)))), Scar ~ Cell,
          value.var = "Presence")
  scar.cell.matrix[is.na(scar.cell.matrix)] <- 0
  # NEW
  all.scar.temp <- aggregate(cs$Cell,
                             by = list(Scar = cs$Scar),
                             length)
  colnames(all.scar.temp)[2] <- "All.scar"
  # OLD
  # all.scar.temp <- aggregate(cs$Cell,
                             # by = list(Scar = cs$Scar,
                                       # Cell.type = cs$Cell.type),
                             # length)
  # colnames(all.scar.temp)[3] <- "All.scar"

  # To calculate how many cells only have a specific scar: these are the cells
  # that have that scar + have only one scar. So if we construct the list of
  # cells that only have one scar and then combine that with which cells have
  # which scars, we should get somewhere?
  
  # NEW
  cell.scar.counts <- aggregate(cs$Scar,
                                by = list(Cell = cs$Cell),
                                length)
  # OLD
  # cell.scar.counts <- aggregate(cs$Scar,
  #                               by = list(Cell = cs$Cell,
  #                                         Cell.type = cs$Cell.type),
  #                               length)

  only.scar <- cs[cs$Cell %in% cell.scar.counts$Cell[cell.scar.counts$x == 1], ]
  if(nrow(only.scar) == 0){
    # NEW
    only.scar.count <- data.frame(Scar = character(),
                                  # Cell.type = character(),
                                  Only.scar = integer())
    # OLD
    # only.scar.count <- data.frame(Scar = character(),
    #                               Cell.type = character(),
    #                               Only.scar = integer())
  }else{
    # NEW
    only.scar.count <- aggregate(only.scar$Cell,
                                 by = list(Scar = only.scar$Scar),
                                 length)
    colnames(only.scar.count)[2] <- "Only.scar"
    # OLD
    # only.scar.count <- aggregate(only.scar$Cell,
    #                              by = list(Scar = only.scar$Scar,
    #                                        Cell.type = only.scar$Cell.type),
    #                              length)
    # colnames(only.scar.count)[3] <- "Only.scar"
  }
  all.scar.temp <- merge(all.scar.temp, only.scar.count, all = T)
  all.scar.temp[is.na(all.scar.temp)] <- 0

  # NEW
  all.scar.temp$Total <- length(unique(cs$Cell))
  # OLD
  # cell.type.cells <- unique(cs[, c("Cell", "Cell.type")])
  # cell.type.counts <- as.data.frame(table(cell.type.cells$Cell.type))
  # colnames(cell.type.counts) <- c("Cell.type", "Total")
  # all.scar.temp <- merge(all.scar.temp, cell.type.counts)
  
  all.scar.temp$Scar.other <-
    all.scar.temp$All.scar - all.scar.temp$Only.scar
  all.scar.temp$Only.other <-
    all.scar.temp$Total - all.scar.temp$All.scar
  all.scar.temp$Total.other <-
    all.scar.temp$Only.other + all.scar.temp$Scar.other
  all.scar.temp$p_A <-
    all.scar.temp$Scar.other/all.scar.temp$Total.other
  all.scar.temp$Log.one.minus.pA <- log(1 - all.scar.temp$p_A)
  all.scar.temp <- all.scar.temp[complete.cases(all.scar.temp), ]

  # possible.top.scars <- 
  #   all.scar.temp$Scar[all.scar.temp$p_A >= max(all.scar.temp$p_A)/10]
  
  # possible.top.scars <- all.scar.temp
  # max.p.celltypes <- aggregate(all.scar.temp$p_A,
  #                              by = list(Cell.type = all.scar.temp$Cell.type),
  #                              max)
  # colnames(max.p.celltypes)[2] <- "Max.p"
  # possible.top.scars <- merge(all.scar.temp, max.p.celltypes)
  # possible.top.scars <- 
  #   possible.top.scars[possible.top.scars$p_A >= possible.top.scars$Max.p/10, ]
  # 
  # scar.lls <- scar.lls[scar.lls$Scar %in% unique(possible.top.scars$Scar), ]
  
  # Construct a scar-scar matrix conn.p that contains the chance of observing a
  # connection x-y for scars x, y, under the assumption that x is the topmost
  # scar. The chance of observing a connection x-y is given by
  # p(O(x<->y)) = 1 - p(!O(x<-y)) = 
  #               1 - Product_{A in celltypes}(1 - p_A(x)^{y_A}).
  # After rewriting log(p(!O(x<->y))) = Sum_{A in celltypes}y_A*log(1 - p_A(x)),
  # we can calculate log(p(!O(x<->y))) as a matrix [log.noconn.p] using matrix 
  # multiplication of y_A [celltype.scar.m] with log(1 - p_A) [lomp.m]. We can
  # then exponentiate log(p(!O(x<->y))) and subtract it from 1 to get
  # p(O(x<->y)) [conn.p]. Since we cannot have connections x<->x, we set the
  # diagonal to 0.
  # This calculation has one complication: if p_A(x) = 1 for a certain A and
  # x, the logarithms produce -Inf. In itself not problematic but the values do 
  # not behave as they should in logspace (p_A(x)=1 should mean that the chance
  # of a connection x<->y is 1 for any y observed in celltype A, and everything
  # else is unaffected). To get the correct outcome, we first set the -Inf in
  # lomp.m to zero to compute everything that should not be affected by 
  # p=1 (i.e. all scars y that do not occur in celltypes A with p_A(x) = 1).
  # Then identify which connections x <-> y are affected and change the p-values
  # for these in conn.p
  
  # NEW
  # Column scar is the scar assumed to be the top scar
  log.noconn.p <- all.scar.temp$All.scar %*% t(all.scar.temp$Log.one.minus.pA)
  colnames(log.noconn.p) <- all.scar.temp$Scar
  rownames(log.noconn.p) <- all.scar.temp$Scar
  conn.p <- 1 - exp(log.noconn.p)
  diag(conn.p) <- 0
  
  # OLD
  # celltype.scar.m <- acast(all.scar.temp[all.scar.temp$p_A != 1, ],
  #                          Scar ~ Cell.type, value.var = "All.scar")
  # celltype.scar.m <- acast(all.scar.temp,
  #                          Scar ~ Cell.type, value.var = "All.scar")
  # celltype.scar.m[is.na(celltype.scar.m)] <- 0
  # lomp.m <- acast(all.scar.temp,
  #                 Cell.type ~ Scar, value.var = "Log.one.minus.pA")
  # lomp.m[is.na(lomp.m)] <- 0

  # Set the lomp.m-value for scars with p=1 to 0 to correctly compute all
  # entries except for the ones that do indeed have other scars in cells of
  # that cell type.
  # lomp.m[lomp.m < 0 & is.infinite(lomp.m)] <- 0
  
  # Compute the connection p-values
  # log.noconn.p <- celltype.scar.m %*% lomp.m
  # log.noconn.p[is.nan(log.noconn.p)] <- 0
  # conn.p <- 1 - exp(log.noconn.p)
  # diag(conn.p) <- 0
  
  # This computation now has incorrect values for scar combinations affected
  # by p_A(x) = 1.
  # sure.shots <- all.scar.temp[all.scar.temp$p_A == 1, c("Scar", "Cell.type")]
  # sure.shot.partners <- merge(cs, sure.shots)
  # colnames(sure.shot.partners)[1] <- "Sure.shot.scar"
  # sure.shot.partners <- merge(sure.shot.partners, cs)
  # sure.shot.partners <- sure.shot.partners[sure.shot.partners$Sure.shot.scar !=
  #                                            sure.shot.partners$Scar, ]
  # for(shot in 1:nrow(sure.shot.partners)){
  #   conn.p[rownames(conn.p) == sure.shot.partners$Scar[shot],
  #          colnames(conn.p) == sure.shot.partners$Sure.shot.scar[shot]] <- 1
  # }
  # END OLD
  
  # Calculate expected degree
  exp.degree <- data.frame(Scar = colnames(conn.p),
                           Expected.degree = colSums(conn.p))
  scar.lls <- merge(scar.lls, exp.degree)
  # Calculate p-value of degree observation
  # x <- scar.lls[2, c("Scar", "Degree", "Expected.degree")]
  scar.lls$Degree.p <-
    apply(scar.lls[, c("Scar", "Degree", "Expected.degree")], 1,
        function(x) {
          scar <- as.character(x[1])
          d_scar <- as.integer(x[2])
          de_scar <- as.numeric(x[3])
          p.conn.scar <- conn.p[, colnames(conn.p) == scar]
          if(d_scar == (nrow(scar.lls) - 1)){
            # If the measured degree is maximal (i.e. the number of other
            # scars) the degree is set to 1 without testing.
            degree.p <- 1
            # d_scar <- d_scar - (d_scar - de_scar)/2
          }else if(d_scar <= de_scar){ # Observation less than expected - lower tail
            degree.p <- ppoisbinom(d_scar, p.conn.scar)
          }else{ # Observation more than expected - higher tail
            degree.p <- ppoisbinom(d_scar, p.conn.scar, lower_tail = F)
          }
          return(degree.p)
        })

  # NEW
  scar.lls <- merge(scar.lls, all.scar.temp[, c("Scar", "Total.other", "p_A")])
  # OLD
  # scar.lls <- merge(scar.lls, all.scar.temp[, c("Scar", "Cell.type", "Total.other", "p_A")])

  scar.lls <- scar.lls[order(-scar.lls$Degree.p, -scar.lls$Degree, 
                             -scar.lls$Scar.count), ]
  
  return(scar.lls)
}

create.ip <- function(readout, scar.probabilities, ip = "Euclidean"){
  scars.with.p <- merge(readout, scar.probabilities[, c("Scar", "p")])
  scars.with.p <- unique(scars.with.p[, c("Scar", "p")])
  scars.with.p$Presence <- 1
  scars.with.p$Logminusp <- -log(scars.with.p$p)
  
  # Create scar metrics
  if(ip == "Euclidean"){
    scar.inner <- acast(scars.with.p, Scar ~ Scar, value.var ="Presence")
    scar.inner[is.na(scar.inner)] <- 0
  }else if (ip %in% c("Sumneglog", "Suminvp")){
    scar.inner <- acast(scars.with.p, Scar ~ Scar, value.var = "Logminusp")
    scar.inner[is.na(scar.inner)] <- 0
  }
  if(ip == "Suminvp"){
    scar.inner.2 <- exp(scar.inner)
    diag(scar.inner) <- diag(scar.inner.2)
  }
  return(scar.inner)
}

create.scar.positions <- function(readout, imputation.method = "none"){
  scars.and.cells <- readout
  scars.and.cells$Presence <- 1
  scars.and.cells.observed <- 
    acast(scars.and.cells, Scar ~ Cell, value.var = "Presence")
  scars.and.cells.observed[is.na(scars.and.cells.observed)] <- 0
  outcell.vector <- data.frame(outcell = rep(0, nrow(scars.and.cells.observed)))
  scars.and.cells.observed.oc <- 
    as.matrix(cbind(scars.and.cells.observed, outcell.vector))
  scars.cells <- scars.and.cells.observed.oc
  
  if(imputation.method %in% c("first-order", "first-correction")){
    # Create scar-cell matrix imputed values
    present <- data.frame(table(readout$Scar))
    colnames(present) <- c("Scar", "Freq")
    # What are the conditional probabilities between scars?
    # "From" is the predicting scar, "To" is the predicted scar.
    # In formula, we are calculating P("To"|"From") = P("To" ^ "From")/P("From")
    scars.coinciding.m <- scars.and.cells.observed %*% t(scars.and.cells.observed)
    scars.coinciding.l <- expand.grid(colnames(scars.coinciding.m), colnames(scars.coinciding.m))
    scars.coinciding.l <- melt(scars.coinciding.m)
    colnames(scars.coinciding.l) <- c("From", "To", "Count") 
    scars.coinciding <- merge(scars.coinciding.l, present, by.x = "From", by.y = "Scar")
    scars.coinciding$CP <- scars.coinciding$Count/scars.coinciding$Freq
    
    # The chance of not observing a scar despite it
    # being present is p(!A) = (1 - p(A|B)) * (1 - p(A|C)) * ... with the product
    # taken over all scars present in the cell.
    # To calculate this, we first observe that for a cell i and a scar A,
    # p^i_A = 1 - exp(log(1 - p^i_A)), and that
    # LNP^i_A := log(1 - p^i_A) = sum_B(S^i_B * log(1 - p(A|B))), which can be easily 
    # calculated as matrix multiplication of the scar detection matrix (S^i_B) with 
    # log(1 minus the conditional scar probabilities).
    # scar.presence.matrix <- t(m.scars)
    scars.coinciding.cpm <- 
      log(1.0000000001 - acast(scars.coinciding, From ~ To, value.var = "CP"))
    # Pseudocount to deal with infinities that otherwise occur.
    LNP <- t(scars.and.cells.observed) %*% scars.coinciding.cpm
    # Note that exp(M) in R is the element-wise exponentiation of the matrix M.
    imputed.scars.pre <- 1 - exp(LNP)
    imputed.scars <- round(imputed.scars.pre, 6)
    scars.and.cells.imputed <- t(imputed.scars)
    scars.and.cells.imputed.oc <- 
      as.matrix(cbind(scars.and.cells.imputed, outcell.vector))
    
    scars.cells <- scars.and.cells.imputed.oc
  }
  
  if(imputation.method == "first-correction"){
    # If scars coincide in at least one cell, their respective conditional
    # probability will be more than zero. For example, if scars A and B co-occur
    # at least once, and a cell has scar A, this cell has a non-zero probability
    # of also having scar B. However, if we also observe scar C in this cell, and
    # scar C never co-occurs with scar B, we'd like to set the probability for B
    # back to zero.
    
    # To do this, we create the matrix of all scars allowed in each cell as
    # follows. If C is the matrix of all scar coincidences, and O the matrix
    # of all scar observations in all cells, then !C is the matrix of all scars
    # that do not coincide even once. If scar A and B do not coincide even once,
    # then the presence of A implies the absence of B. This absence is encoded
    # in !C * O, and the matrix of allowed scars in each cell therefore is
    # given by Allow = !(!C * O).
    
    # Finally, we do an elementwise multiplication of the imputation values
    # with the matrix of allowed scars to get the imputation values of allowed
    # scars.
    neg.scars.coinciding <- !(scars.coinciding.m)
    scars.excluded <- neg.scars.coinciding %*% scars.and.cells.observed.oc
    scars.allowed <- !(scars.excluded)
    scars.cells <- scars.cells * scars.allowed
  }
  
  if(imputation.method == "full"){
    full.check <- t(scars.cells) %*% scars.cells
    # c <- 1
    # f.imputed <- 0
    for(c in 1:(ncol(scars.cells) - 1)){
      # print(c)
      cells.rel.subset.1 <- rownames(full.check)[full.check[, c] == max(full.check[, c])]
      cells.rel.subset <- setdiff(cells.rel.subset.1, colnames(scars.cells)[c])
      
      if(length(cells.rel.subset) == 0){next}
      # f.imputed <- f.imputed + 1
      scars.cells.rel.subset <- 
        scars.cells[, colnames(scars.cells) %in% cells.rel.subset]
      if(length(cells.rel.subset) == 1){
        scars.cells[, c] <- scars.cells.rel.subset
      }else{
        scars.cells[, c] <- rowMeans(scars.cells.rel.subset)
      }
    }
  }
  
  return(scars.cells)
}

cull.readout <- function(readout, freq.prod){
  scar.freqs <- data.frame(table(readout$Scar))
  colnames(scar.freqs) <- c("Scar", "Cells")
  scar.freqs$Frequency <- scar.freqs$Cells/length(unique(readout$Cell))
  readout.plus <- merge(readout, scar.freqs[, c("Scar", "Frequency")])
  cell.scar.freqs <-
    aggregate(readout.plus$Frequency,
              by = list(Cell = readout.plus$Cell),
              prod)
  colnames(cell.scar.freqs)[2] <- "Freq.prod"
  cells.keep <- cell.scar.freqs$Cell[cell.scar.freqs$Freq.prod < freq.prod]
  readout.final <- readout[readout$Cell %in% cells.keep, ]
  
  return(readout.final)
}

decimal <- function(x){
  output <- 0
  power <- 0
  while(x != ""){
    output <- output + as.integer(substr(x, nchar(x), nchar(x))) * 2 ^ power
    x <- substr(x, 1, nchar(x) - 1)
    power <- power + 1
  }
  
  return(output)
}

determine.main <- function(aggregates, current.branches, main.parent, 
                           branch.size.ratio){
  # Recursive function: calculate what's main and off-main for a given dataframe of
  # branches (these should be sister clades), determine main and off-main for
  # children, return dataframe with main and off-main marked.
  # Create output dataframe
  output <- current.branches
  output$Main <- F
  # main.parent F -> Set all branches off main (happens by default)
  if(main.parent){
    # main.parent T -> Test if there are any branches off-main
    if(nrow(output) > 1){
      output$Max.freq <- max(output$Freq)
      output$Freq.ratio <- output$Freq/output$Max.freq
      output$Main <- output$Freq.ratio >= branch.size.ratio
      
       output <- output[, c("Node", "Node.depth", "Freq", "Main")]
    }else{
      output$Main <- T
    }
  }
  # For each branch: 
  for(i in 1:nrow(output)){
    #   determine children
    parent <- output$Node[i]
    children <- 
      aggregates[grepl(paste(parent, "_", sep = ""), aggregates$Node) &
                   aggregates$Node.depth == (output$Node.depth[i] + 1), ]
    if(nrow(children) > 0){
      #   if children -> call function for children, add children to output
      children.main <-
        determine.main(aggregates, current.branches = children, 
                       main.parent = output$Main[i], 
                       branch.size.ratio = branch.size.ratio)
      output <- rbind(output, children.main)
    }
  }
  # return output
  return(output)
}

expand.tip <- function(x, y, attach.at){
  # Making a new tree z that is basically tree x with one tip replaced by a
  # second tree y.
  # The new tree will have Nnode(z) = Nnode(x) + Nnode(y).
  # We will leave out the node labels.
  attach.tip <- which(x$tip.label == attach.at)
  
  z.edges <- x$edge
  aux.node.number <- max(x$edge) + 1
  z.edges[z.edges == attach.tip] <- aux.node.number
  z.edges[z.edges > attach.tip] <- z.edges[z.edges > attach.tip] - 1
  z.edges <- z.edges + length(y$tip.label)
  y.edges <- y$edge
  y.edges[y.edges > length(y$tip.label)] <-
    y.edges[y.edges > length(y$tip.label)] + x$Nnode + length(x$tip.label) - 1
  
  z.edges <- rbind(z.edges, y.edges)
  z <- list(edge = z.edges,
            tip.label = c(y$tip.label, x$tip.label[-attach.tip]),
            Nnode = x$Nnode + y$Nnode)
  
  if(!is.null(x$edge.length) & !is.null(y$edge.length)){
    z.edge.lengths <- c(x$edge.length, y$edge.length)
    z$edge.length <- z.edge.lengths
  }
  
  class(z) <- "phylo"
  
  return(z)
}

# Getting scar indices
f <- function (i, j, n) {
  ifelse((i > j) & (j <= n), (j - 1) * (2 * n - 2 - j) / 2 + (i - 1), NA_real_)
}

generate_tree = function(df){
  # columns.include <- c("Parent", "Child", "Scar.acquisition")
  # if(!is.null(fill.col)){
  #   columns.include <- c(columns.include, fill.col)
  # }
  # if(!is.null(size.col)){ columns.include <- c(columns.include, size.col)}
  for(i in 1:nrow(df)){
    parent = paste0('nd', as.character(df$Parent[i]))
    child = paste0('nd', as.character(df$Child[i]))
    scar <- df$Scar.acquisition[i]
    cell.type <- df$Cell.type[i]

    if(!exists(child)){
      eval_txt = sprintf('%s <<- Node$new("%s", name="%s", scar = "%s", Cell.type = "%s")',
                         child, child, child, scar, cell.type)
      eval(parse(text=eval_txt))
      if("fill" %in% colnames(df)){
        eval_txt <- paste(child, "$fill <- \"", df$fill[i], "\"", sep = "")
        eval(parse(text = eval_txt))
      }
      if("size" %in% colnames(df)){
        eval_txt <- paste(child, "$size <- ", df$size[i], sep = "")
        eval(parse(text = eval_txt))
      }
      if("Main" %in% colnames(df)){
        eval_txt <- paste(child, "$Main <- ", df$Main[i], sep = "")
        eval(parse(text = eval_txt))
      }
    }
    if(exists(parent)){
      add_txt = sprintf('%s$AddChildNode(%s)', parent, child)
      eval(parse(text=add_txt))
    }
  }

  return_tree <- eval(parse(text=sprintf('%s$root', ls(envir=globalenv(), pattern='^nd')[1])))
  rm(list=ls(envir=globalenv(), pattern='^nd'), envir=globalenv())

  return(return_tree)
}


get.dist.index <- function(k, n){
  # Get indices for the dth element of a distance object of distances between
  # n elements.
  # Starting position for each column
  ptr_all_cols <- f(2:n, 1:(n - 1), n)
  # Maximum valid `k`
  k_max <- n * (n - 1) / 2
  
  if (k > k_max) return(c(i = NA_real_, j = NA_real_))
  j <- sum(ptr_all_cols <= k)  ## get column index j
  i <- k - ptr_all_cols[j] + j + 1  ## get row index i
  return(c(i = i, j = j))
}

filter.data <- function(X.log, expression.threshold = 3, gene.threshold = 300,
                        filter.genes = F){
  genes.detected <-
    data.frame(Cell = colnames(X.log),
               Genes = 
                 apply(X.log, 2,
                       # function(x) sum(x >= log(expression.threshold))))
                       function(x) sum(exp(x) >= (expression.threshold + 1))))
  genes.in.cells <-
    data.frame(Gene = rownames(X.log),
               Cells =
                 apply(X.log, 1, 
                       # function(x) sum(x >= log(expression.threshold))))
                       function(x) sum(exp(x) >= (expression.threshold + 1))))
  genes.in.cells$Keep <- T
  if(filter.genes){
    genes.in.cells$Keep <- (genes.in.cells$Cells > 0)
  }
  X.out <- X.log[genes.in.cells$Keep, genes.detected$Genes >= gene.threshold]
  
  return(X.out)
}

filter.data.lin <- function(X, expression.threshold = 3, gene.threshold = 300,
                        filter.genes = F){
  genes.detected <-
    data.frame(Cell = colnames(X),
               Genes = 
                 apply(X, 2,
                       function(x) sum(x >= expression.threshold)))
  genes.in.cells <-
    data.frame(Gene = rownames(X),
               Cells =
                 apply(X, 1, 
                       function(x) sum(x >= expression.threshold)))
  genes.in.cells$Keep <- T
  if(filter.genes){
    genes.in.cells$Keep <- (genes.in.cells$Cells > 0)
  }
  X.out <- X[genes.in.cells$Keep, genes.detected$Genes >= gene.threshold]
  
  return(X.out)
}

find.all.p.links <- function(scars.w.p){
  if(nrow(scars.w.p) == 0){
    return(data.frame(Cell.1 = character(),
                      Cell.2 = character(),
                      p = numeric()))
  }
  # scars.w.p$psq <- sqrt(scars.w.p$p)
  scars.w.p$sqrtlogminusp <- sqrt(-log(scars.w.p$p))
  
  cell.scar <- dcast(scars.w.p, Scar ~ Cell, value.var = "sqrtlogminusp")
  cell.scar[is.na(cell.scar)] <- 0
  rownames(cell.scar) <- cell.scar$Scar
  cell.scar <- as.matrix(cell.scar[, -1])
  
  # cell.links <- t(cell.scar) %*% cell.scar
  cell.links.1 <- t(cell.scar) %*% cell.scar
  cell.links <- 1/exp(cell.links.1)
  
  # Links between a cell and itself aren't well-defined. Our convention is to 
  # set them to zero.
  diag(cell.links) <- 0
  
  connections <- data.frame(melt(cell.links,
                                 varnames = c("Cell.1", "Cell.2"),
                                 value.name = "Strength"))
  connections <- unique(connections[!(connections$Strength %in% c(0, 1)),])
  connections$Cell.1 <- as.character(connections$Cell.1)
  connections$Cell.2 <- as.character(connections$Cell.2)

  return(connections)
}

find.scars.to.remove <- function(starting.scar, tree.summary){
  # Find all scars preceding a given [starting.scar], including that scar
  # itself, in a phylogenetic dataframe [tree.summary].
  if(starting.scar == "Root"){scars.to.remove = character()
  }else{
    scars.to.remove <- starting.scar
    preceding.scars <- 
      unique(c(scars.to.remove, 
               as.character(tree.summary$Node.1[tree.summary$Node.2 %in% 
                                                  scars.to.remove])))
    while(!("Root" %in% preceding.scars)){
      scars.to.remove <- preceding.scars
      preceding.scars <- 
        unique(c(scars.to.remove, 
                 as.character(tree.summary$Node.1[tree.summary$Node.2 %in% 
                                                    scars.to.remove])))
      # tree.summary$Node.1[tree.summary$Node.2 %in% scars.to.remove]
    }
  }
  
  return(scars.to.remove)
}

find.weighted.links <- function(fscars.r, method = "log_inv"){
  # Given a data frame with cells and scars, return a dataframe where
  # every line is a unique connection (made by one or more scars) between
  # two cells, together with its weight.
  if(nrow(fscars.r) == 0){
    return(data.frame(Cell.1 = character(),
                      Cell.2 = character()))
  }
  if(method == "log_inv"){
    fscars.r$Strength <- sqrt(log(1/(fscars.r$p)))
  }else if (method == "linear"){
    fscars.r$Strength <- sqrt(fscars.r$p)
  }
  cell.scar <- dcast(fscars.r, Scar ~ Cell, value.var = "Strength")
  cell.scar[is.na(cell.scar)] <- 0

  rownames(cell.scar) <- cell.scar$Scar
  cell.scar <- as.matrix(cell.scar[, -1])
  cell.links <- t(cell.scar) %*% cell.scar
  
  # Links between a cell and itself aren't well-defined. Our convention is to 
  # set them to zero. To avoid double counting we also set the upper triangle
  # of the links to zero.
  diag(cell.links) <- 0

  cell.links[upper.tri(cell.links, diag = T)] <- 0
  connections <- data.frame(melt(cell.links,
                                 varnames = c("Cell.1", "Cell.2"),
                                 value.name = "Strength"))
  # Note that since we set the diagonal and upper triangle to zero, removing 
  # all links with strength 0 automatically removes self-connections and double
  # counting.
  connections <- unique(connections[connections$Strength != 0,])
  connections$Cell.1 <- as.character(connections$Cell.1)
  connections$Cell.2 <- as.character(connections$Cell.2)
  
  return(connections)
}

find.links <- function(fscars.r){
  # Given a data frame with cells and scars, return a dataframe where
  # every line is a unique connection (made by one or more scars) between
  # two cells. 
  if(nrow(fscars.r) == 0){
    return(data.frame(Cell.1 = character(),
                      Cell.2 = character()))
  }
  cell.scar <- dcast(fscars.r, Scar ~ Cell, fun.aggregate = length,
                     value.var = "Scar")
  
  cell.scar <- cell.scar[complete.cases(cell.scar), ]
  rownames(cell.scar) <- cell.scar$Scar
  cell.scar <- as.matrix(cell.scar[, -1])
  cell.links <- t(cell.scar) %*% cell.scar
  
  # So far, this has been matrix multiplication, but links are a binary entity
  # and links between a cell and itself aren't possible. So we remove the 
  # diagonal and binarize the rest.
  diag(cell.links) <- 0
  cell.links[cell.links > 0] <- 1
  
  cell.links[upper.tri(cell.links)] <- 0
  connections <- data.frame(melt(cell.links,
                                 varnames = c("Cell.1", "Cell.2"),
                                 value.name = "Link"))
  colnames(connections)[3] <- "Link"
  connections <- unique(connections[connections$Link == 1, -3])
  connections$Cell.1 <- as.character(connections$Cell.1)
  connections$Cell.2 <- as.character(connections$Cell.2)
  
  return(connections)
}

get.readout <- function(scar.cells.final, cells.sampled, 
                        integration.sites, doublet.rate = 0){
  sites <- nrow(integration.sites)
  
  # cells.included = singlets + 2 * doublets = cells.sampled + number.doublets
  number.doublets <- rbinom(1, cells.sampled, doublet.rate)
  cells.included <- cells.sampled + number.doublets
  
  cells.found <-
    data.frame(Cell = sample(scar.cells.final$Cell, cells.included))
  cells.found <- merge(cells.found, scar.cells.final[, c("Cell", "Cell.type", 
                                                         "Detection.rate")])
  cells.found$Scars.found <-
    rbinom(n = cells.included, size = sites, prob = cells.found$Detection.rate)
  full.readout <- scar.cells.final[scar.cells.final$Cell %in% cells.found$Cell, ]
  long.drop.out.readout <-
    data.frame(Cell = rep(cells.found$Cell, cells.found$Scars.found),
               Scar = NA)
  long.drop.out.readout <- merge(long.drop.out.readout, 
                                 scar.cells.final[, c("Cell", "Cell.type")])
                                 
  
  # First determine which integration sites are found in a cell, then
  # determine whether we actually see those integration sites (using the
  # integration site detection rates).
  readout.row <- 1
  c.row <- 1 # For debugging
  for(c.row in 1:cells.included){
    c.scar.size <- cells.found$Scars.found[c.row]
    if(c.scar.size == 0){next}
    cell <- cells.found$Cell[c.row]
    scars.this.cell <- 
      melt(sample(full.readout[full.readout$Cell == cell, 2 + 1:sites], 
                  c.scar.size), id.vars = NULL)
    colnames(scars.this.cell) <- c("Site", "Scar")
    scars.this.cell <- merge(scars.this.cell, integration.sites)
    scars.this.cell$Detection.roll <- runif(nrow(scars.this.cell))
    scars.this.cell$Detected <- scars.this.cell$Detection.rate >= scars.this.cell$Detection.roll
    scars.this.cell$Scar <- ifelse(scars.this.cell$Detected, scars.this.cell$Scar, NA)
    long.drop.out.readout$Scar[readout.row:(readout.row + c.scar.size - 1)] <-
      scars.this.cell$Scar
    readout.row <- readout.row + c.scar.size
  }
  
  readout.pre.doublet <- long.drop.out.readout[complete.cases(long.drop.out.readout), ]
  
  readout.pre.doublet <- unique(readout.pre.doublet)
  
  if(doublet.rate > 0){
    doublets <- 
      data.frame(
        matrix(cells.found$Cell[sample.int(nrow(cells.found), 
                                           2 * number.doublets)], ncol = 2))
    doublets$Cell <- paste(doublets$X1, doublets$X2, sep = ";")
    doublet.scars <- 
      rbind(merge(doublets, readout.pre.doublet[, c("Cell", "Scar")], by.x = "X1", 
                  by.y = "Cell")[, c("Cell", "Scar")],
            merge(doublets, readout.pre.doublet[, c("Cell", "Scar")], by.x = "X2", 
                  by.y = "Cell")[, c("Cell", "Scar")])
    doublet.scars$Cell.type <- "Doublet"
    # The total number of doublets can be lower than the number of doublets
    # we started with since not all cells have a scar readout.
  readout.no.doublet <- 
    readout.pre.doublet[!(readout.pre.doublet$Cell %in% 
                            c(as.character(doublets$X1), 
                              as.character(doublets$X2))), ]
    readout <- rbind(readout.no.doublet, doublet.scars)
  }else{
    readout <- readout.pre.doublet
  }
  readout <- unique(readout)
  
  return(readout)
}

get.scar.names <- function(fscars, wt.remove = T, return.seq.and.loc = F,
                           wt.seq = wildtype.seq){
  # Return a dataframe with scar names and remove wildtype if wt.remove is True.
  present.scars <- unique(fscars[, c("CIGAR", "Sequence", "Location")])
  present.scars <- merge(present.scars, scar.probabilities[, c("Sequence", "p")])
  
  
  present.scars$Scar <- 
    paste(1:nrow(present.scars), present.scars$CIGAR, sep = ":")
  fscars <-
    merge(fscars, present.scars)
  if(return.seq.and.loc){
    fscars.r <- unique(fscars[, c("Cell", "Scar", "Location", "Sequence")])
  }else{
    fscars.r <- unique(fscars[, c("Cell", "Scar")])
  }
  
  if(wt.remove){
    wildtype <- present.scars$Scar[present.scars$Sequence == wt.seq]
    fscars.r <- fscars.r[fscars.r$Scar != wildtype, ]
  }
  
  return(fscars.r)
}

get.connections <- function(scars, scars.present.probabilities, wildtype.clip, tsne.coord){
  scar.adj <- 
    calculate.adjacency(scars, scars.present.probabilities, wildtype.clip, wt.remove = F)
  
  # NOTE: the connections thus calculated are unique and the assignment of a
  # cell to cell.1 or cell.2 is meaningless.
  scar.adj[lower.tri(scar.adj)] <- 1
  connections <- melt(as.matrix(scar.adj))
  colnames(connections) <- c("Cell.1", "Cell.2", "p")
  connections <- connections[!(connections$p %in% c(0, 1)), ]
  connections <- merge(connections, tsne.coord[, 1:3],
                       by.x = "Cell.1", by.y = "Cell")
  colnames(connections)[4:5] <- 
    paste("Cell.1", colnames(connections)[4:5], sep = ".")
  connections <- merge(connections, tsne.coord[, 1:3],
                       by.x = "Cell.2", by.y = "Cell")
  colnames(connections)[6:7] <-
    paste("Cell.2", colnames(connections)[6:7], sep = ".")
  connections <- 
    merge(connections,
          tsne.coord[, c("Cell", "Cell.type")],
          by.x = "Cell.1", by.y = "Cell")
  colnames(connections)[8] <- "Cell.type.1"
  connections <-
    merge(connections, tsne.coord[, c("Cell", "Cell.type")], 
          by.x = "Cell.2", by.y = "Cell")
  colnames(connections)[9] <- "Cell.type.2"
  
  return(connections)
}

graph.and.decompose <- function(cs){
  # Create scar graph based on scars coinciding in at least one cell
  # and decompose into connected components. Output the
  # list of connected components.
  if(length(unique(cs$Scar)) == 1){
    only.once.connections <- data.frame(Scar.A = character(),
                                        Scar.B = character())
  }else{
    scar.connections <- connections.for.graph(cs)
    only.once.connections <- data.frame(t(combn(unique(cs$Scar), 2)))
    colnames(only.once.connections) <- c("Scar.A", "Scar.B")
    only.once.connections <- 
      merge(only.once.connections, 
            scar.connections[scar.connections$x_AB > 0, ])
    
    # Debugging code START
    # ooc.2 <- only.once.connections
    # ooc.2$AB.ratio.obs <- ooc.2$x_AB/length(unique(cs$Cell))
    # ooc.2$AB.ratio.exp <- ooc.2$x_A * ooc.2$x_B/(length(unique(cs$Cell)))^2
    # ooc.2$AB.obs.exp <- ooc.2$AB.ratio.obs/ooc.2$AB.ratio.exp
    # ooc.cutoff <- ooc.2[ooc.2$AB.obs.exp > 1, ]
    # ooc.cutoff.graph <-
    #   graph_from_data_frame(ooc.cutoff[, c("Scar.A", "Scar.B")],
    #                         directed = F, vertices = union(ooc.cutoff$Scar.A, ooc.cutoff$Scar.B))
    #plot(ooc.cutoff.graph)
    # Debugging code END


  }
  current.full.scar.graph <- 
    graph_from_data_frame(only.once.connections[, c("Scar.A", "Scar.B")],
                          directed = F, vertices = unique(cs$Scar))
  
  decomposed.graph <- decompose(current.full.scar.graph) 
  
  return(decomposed.graph)
}

initialize.branches <- function(current.cs, scar.remove = "Root", 
                                size.ratio = 0.25){
  # Return a dataframe of all branches present in a cell/scar dataset.
  # Mark branches that are less than size.ratio of the next bigger branch
  # as 'not main' (underlying thought is that we missed a scar connection
  # here and therefore won't be able to place that branch correctly).
  
  # Create and decompose scar graph
  decomposed.graph <- graph.and.decompose(current.cs)
  current.comps <- length(decomposed.graph)
  
  # Initialize dataframe of branches
  tree.summary <- data.frame(Node.1 = rep(scar.remove, current.comps),
                             Node.2 = rep(NA, current.comps),
                             Component = 1:current.comps,
                             Depth = rep(0, current.comps),
                             Size = rep(NA, current.comps),
                             Main = rep(F, current.comps),
                             stringsAsFactors = F)
  
  # Determine branch sizes
  for(tree.row in 1:nrow(tree.summary)){
    cc <- tree.summary$Component[tree.row]
    current.graph <- decomposed.graph[[cc]]
    current.cs.component <- current.cs[current.cs$Scar %in% V(current.graph)$name, ]
    tree.summary$Size[cc] <-   length(unique(current.cs.component$Cell))
  }
  tree.summary <- tree.summary[order(-tree.summary$Size), ]
  
  # Mark any branches more than size.ratio smaller than its predecessor as
  # non-main.
  tree.summary$Main[1] <- T
  if(nrow(tree.summary) > 1){
    tree.summary$Previous.size <- c(tree.summary$Size[1], tree.summary$Size[-nrow(tree.summary)])
    tree.summary$Size.ratio <- tree.summary$Size/tree.summary$Previous.size
    for(tree.row in 2:nrow(tree.summary)){
      if(tree.summary$Size.ratio[tree.row] < size.ratio){
        break
      }else{
        tree.summary$Main[tree.row] <- T
      }
    }
  }
  
  return(tree.summary[, c("Node.1", "Node.2", "Component", "Depth", "Size", "Main")])
}

link.fisher.test <- function(cells, connections){
  # Add strata to the connections data
  stratum <- colnames(cells)[2]
  connections <- merge(connections, cells, by.x = "Cell.1", by.y = "Cell")
  colnames(connections)[3] <- paste(stratum, ".1", sep = "")
  connections <- merge(connections, cells, by.x = "Cell.2", by.y = "Cell")
  colnames(connections)[4] <- paste(stratum, ".2", sep = "")
  connections <- connections[, c(2, 1, 3, 4)]
  
  # Prepare the output data
  strata <- unique(cells[ ,2])
  fisher.outcomes <- data.frame(matrix(ncol = 4, nrow = length(strata) ** 2))
  colnames(fisher.outcomes) <- 
    c(paste(stratum, ".1", sep = ""),
      paste(stratum, ".2", sep = ""),
      "Odds.ratio", "p")
  
  i <- 1
  for(stratum.1 in strata){
    for(stratum.2 in strata){
      # Prepare names for this particular test combination
      fisher.outcomes[i, 1] <- stratum.1
      fisher.outcomes[i, 2] <- stratum.2
      test.strata <- c(stratum.1, stratum.2)
      
      # Calculate the true and false positives directly from the connections
      connections$Correct <- 
        (connections[, 3] == stratum.1 & connections[, 4] == stratum.2) |
        (connections[, 3] == stratum.2 & connections[, 4] == stratum.1)
      connections$Incorrect <- 
        ((connections[, 3] == stratum.1 & connections[, 4] != stratum.2) |
           (connections[, 3] == stratum.2 & connections[, 4] != stratum.1) |
           (connections[, 3] != stratum.1 & connections[, 4] == stratum.2) |
           (connections[, 3] != stratum.2 & connections[, 4] == stratum.1))
      true.pos <- sum(connections$Correct)
      false.pos <- sum(connections$Incorrect)
      
      # Calculate the false negatives by calculating all possible true
      # connections
      cells.stratum.1 <- sum(cells[, 2] == stratum.1)
      cells.stratum.2 <- sum(cells[, 2] == stratum.2)
      all.possible.correct.connections <- 
        ifelse(stratum.1 == stratum.2, 
               0.5 * cells.stratum.1 * (cells.stratum.1 - 1),
               cells.stratum.1 * cells.stratum.2)
      false.neg <- all.possible.correct.connections - true.pos
      
      # Calculate the true negatives by first calculating all possible
      # connections, and then calculating all possible false connections.
      all.connections.all.cells <-
        0.5 * (nrow(cells)) * (nrow(cells) - 1)
      cells.outside <- nrow(cells[!(cells[, 2] %in% test.strata), ])
      all.connections.outside <-
        0.5 * cells.outside * (cells.outside - 1)
      all.possible.relevant.connections <-
        all.connections.all.cells - all.connections.outside
      all.possible.incorrect.connections <- 
        all.possible.relevant.connections - all.possible.correct.connections
      true.neg <- all.possible.incorrect.connections - false.pos
      
      # Do Fisher's exact test
      fisher.matrix <- matrix(c(true.pos, false.neg, false.pos, true.neg), 
                              ncol = 2, nrow = 2)
      fisher.test.outcome <- fisher.test(fisher.matrix)
      
      # Add results to the output data
      fisher.outcomes$Odds.ratio[i] <- fisher.test.outcome$estimate
      fisher.outcomes$p[i] <- fisher.test.outcome$p.value
      i <- i + 1
    }
  }
  
  return(fisher.outcomes)
}

make.edgelist <- function(tree.summary, node.counts, cells, main.min = 5,
                          off.main.min = 10){
  # Create an edgelist to generate a tree. All cells are placed in the
  # edgelist.
  
  # Edges without cells first.
  tree.plot <- tree.summary[, c("Node", "Node.2")]
  colnames(tree.plot) <- c("Child", "Scar.acquisition")
  tree.plot$Parent <-
    sapply(tree.plot$Child,
           function(x){
             y <- unlist(strsplit(x, "_"))
             z <- paste(y[-length(y)], collapse = "_")
           }
    )
  root.scars <-
    unique(tree.summary$Node.1[grepl("Root", tree.summary$Node.1)])
  if(root.scars != "Root"){
    root.add <- data.frame(Parent = "Root",
                           Child = 0,
                           Scar.acquisition = root.scars,
                           stringsAsFactors = F)
    root.add$Scar.acquisition <- 
      sapply(root.add$Scar.acquisition,
             function(x) paste(unlist(strsplit(x, ","))[-1], collapse = ","))
    tree.plot <- rbind(root.add, tree.plot)
  }
  tree.plot$Cell.type <- "NA"
  tree.plot$fill <- "black"
  tree.plot$size <- 1
  
  tree.plot <- merge(tree.plot, node.counts[, c("Node", "Main")],
                     by.x = "Child", by.y = "Node")
  if(sum(tree.plot$Parent == "0") > 1){
    tree.plot <-
      rbind(data.frame(Child = 0, Scar.acquisition = "", Parent = "Root", Cell.type = "NA",
                       fill = "black", size = 1, Main = T),
            tree.plot)
    tree.plot$Parent <- as.character(tree.plot$Parent)
  }
  tree.plot$Scar.acquisition <- as.character(tree.plot$Scar.acquisition)
  
  # Add cells
  if("Cell.type" %in% names(cells)){
    cells.add <- cells[, c("Node", "Cell", "Cell.type")]
  }else{
    cells.add <- 
      cells[, c("Node", "Cell")]
    cell.add$Cell.type <- "Cell"
  }
  cells.add$Scar.acquisition <- ""
  colnames(cells.add)[1:2] <- c("Parent", "Child")
  cells.add$fill <- "lightgrey"
  cells.add$size <- 0.5
  cells.add$Child <- 
    sapply(cells.add$Child,
           function(x){
             x <- as.character(x)
             if(grepl(";", x)){
               return(paste(unlist(strsplit(x, ";")), collapse = "d"))
             }else{
               return(x)
             }
           }
    )
  cells.add$Main <- T
  
  # Set thresholds to plot nodes - try plotting nodes >= 5 for main, > 10 for off-main
  # Note - has to hold for the full information plot as well, there we need to collapse
  # the scars as well (node sizes remain the same).
  node.converter <- node.counts
  
  node.converter$To <-
    sapply(node.converter$Node,
           function(x){
             y <- unlist(strsplit(x, "_"))
             return(paste(y[-length(y)], collapse = "_"))
           }
    )
  node.converter$Keep <-
    ifelse(node.converter$Main, 
           node.converter$Freq >= main.min,
           node.converter$Freq >= off.main.min)
  # See if any children of keep = F nodes still have keep = T.
  keep.f.pattern <- paste(node.converter$Node[!node.converter$Keep], collapse = "_|")
  node.converter$Keep[node.converter$Keep] <- !grepl(keep.f.pattern, node.converter$Node[node.converter$Keep])

  node.converter$To <- ifelse(node.converter$Keep, node.converter$Node, node.converter$To)
  
  node.converter.2 <- node.converter[, c("Node", "To")]
  colnames(node.converter.2) <- c("To", "To.2")
  
  while(sum(node.converter$To %in% c("0", node.converter$Node[node.converter$Keep])) < 
        nrow(node.converter)){
  node.converter <- merge(node.converter, node.converter.2, all.x = T)
  node.converter$To.2[is.na(node.converter$To.2)] <- 
    node.converter$To[is.na(node.converter$To.2)]
  node.converter <- node.converter[, -1]
  colnames(node.converter)[which(colnames(node.converter) == "To.2")] <- "To"
  }
  
  if(sum(node.converter$To %in% c("0", node.converter$Node[node.converter$Keep])) < 
     nrow(node.converter)){
    stop("Was not able to collapse small nodes")
  }
  
  # Remove nodes that have keep = F from nodelist.
  tree.plot <- tree.plot[tree.plot$Child %in% c("0", node.converter$Node[node.converter$Keep]), ]
  
  # Shift cells
  cells.add <- merge(cells.add, node.converter[, c("Node", "To")],
                     by.x = "Parent", by.y = "Node")
  cells.add$Parent <- cells.add$To
  cells.add <- cells.add[, -which(colnames(cells.add) == "To")]
  
  # Output
  tree.plot.cells <- rbind(tree.plot, cells.add)
  return(tree.plot.cells)
}

normalize.data <- function(X.log.f, method, input){
  if(method == "median" && input == "log"){
    X.f <- exp(X.log.f)
    med.expression <- median(colSums(X.f))
    X.out.lin <- med.expression * t(t(X.f)/colSums(X.f))
    X.out <- data.frame(log(X.out.lin))
  }
  
  if(method == "median" && input == "lin"){
    X.f <- X.log.f
    med.expression <- median(colSums(X.f))
    X.out.lin <- data.frame(med.expression * t(t(X.f)/colSums(X.f)))
    X.out <- X.out.lin
  }
  
  return(X.out)
}

overseq <- function(x, y, name, plot = T, xlimit = 5){
  os <- as.matrix(x/y)
  os.melt <- melt(os)
  os.melt$Log <- log2(os.melt$value)
  os.melt <- os.melt[is.finite(os.melt$Log), ] # removing all infinities
  # created by a gene having no UMIs in a cell.
  if(plot){
    suppressWarnings(
      print(
        ggplot(os.melt, aes(x = Log)) + geom_histogram(binwidth = 0.1) +
          labs(title = paste("Oversequencing", name),
               x = "Read counts / UMI counts (log2)") +
          scale_x_continuous(breaks = 0:xlimit, limits = c(-0.1, xlimit))
      )  
    )
  }
  
  ovs <- mean(os.melt$value)
  return(ovs)  
}

plot.scars <- function(scars, strata, scatter = T){
  scars.s1 <- scars[scars[, 3] == strata[1], ]
  scars.s1.pl <-
    aggregate(scars.s1$Scar,
              by = list(Scar = scars.s1$Scar),
              length)
  colnames(scars.s1.pl)[2] <- "S1"
  
  scars.s2 <- scars[scars[, 3] == strata[2], ]
  scars.s2.pl <-
    aggregate(scars.s2$Scar,
              by = list(Scar = scars.s2$Scar),
              length)
  colnames(scars.s2.pl)[2] <- "S2"
  
  scars.plot <- merge(scars.s1.pl, scars.s2.pl, all = T)
  scars.plot[is.na(scars.plot)] <- 0
  
  scars.plot$Total <- scars.plot$S1 + scars.plot$S2
  scars.plot$Link.1 <- choose(scars.plot$S1, 2)
  scars.plot$Link.2 <- choose(scars.plot$S2, 2)
  scars.plot$Link.12 <- scars.plot$S1 * scars.plot$S2
  
  scars.plot$Imp <- (scars.plot$S1 + scars.plot$S2 > 5)
  # scars.plot.agg <- melt(scars.plot[scars.plot$Imp, c(1, 4, 5, 6)])
  scars.plot.ni <- scars.plot[!scars.plot$Imp, ]
  scars.plot.ni.a <-
    data.frame(Scar = "Other",
               Link.1 = sum(scars.plot.ni$Link.1),
               Link.2 = sum(scars.plot.ni$Link.2),
               Link.12 = sum(scars.plot.ni$Link.12))
  scars.plot.agg <-
    rbind(melt(scars.plot[scars.plot$Imp, c(1, 5:7)], id.vars = "Scar"),
          melt(scars.plot.ni.a, id.vars = "Scar"))
  
  if(scatter){
    print(
      ggplot(scars.plot) +
        geom_point(aes(x = S1, y = S2)) +
        labs(x = strata[1], y = strata[2])
    )
  }
  print(
    ggplot(scars.plot.agg) +
      geom_bar(aes(x = variable, fill = Scar, y = value), stat = "identity") +
      scale_y_continuous(limits = c(0, 2000))
  )
  
  return(scars.plot)
}

plot.gene <- function(mapping, expression, gene, plot.title, point.size = 2){
  gene.expression <-
    data.frame(Cell = colnames(expression),
               Expression = as.vector(t(expression[
                 rownames(expression) == gene, ])))
  plot.frame <- merge(mapping, gene.expression)
  print(
    ggplot(plot.frame) +
      geom_point(aes(x = X1, y = X2, color = Expression), size = point.size) +
      labs(title = plot.title, x = "", y = "", color = "") +
      scale_color_gradientn(colours = c("steelblue", "white", "yellow", "red"), 
                            values = rescale(c(-1, 0, 1, 3))) +
      theme(axis.ticks = element_blank(), 
            axis.text = element_blank(),
            panel.grid = element_blank())
  )
}

plot.scar.dendro.heatmap <- function(celltree, positions, plot.top, 
                                     scar.plot.size = -1, lC = ""){
  all.scar.order <- data.frame(Mean = rowMeans(positions))
  all.scar.order$Scar <- rownames(all.scar.order)
  all.scar.order <- all.scar.order[order(-all.scar.order$Mean), ]
  all.scar.order$Cummeanfraction <- cumsum(all.scar.order$Mean)/sum(all.scar.order$Mean)
  scars.keep <- all.scar.order$Scar[all.scar.order$Cummeanfraction <= plot.top]
  
  keep.positions <- positions[rownames(positions) %in% scars.keep, ]

  scar.hc <- as.dendrogram(hclust(dist(keep.positions), method = "centroid"))
  scar.pred.order <- rowMeans(keep.positions)
  scar.hc.ro <- reorder(scar.hc, sort(-scar.pred.order))
  
  if(scar.plot.size <= 0){
    heatmap.2(keep.positions, trace = "none",
              breaks = 101, symbreaks = F, col = colors.hm, key = F,
              dendrogram = "column", Colv = celltree, Rowv = scar.hc.ro,
              labRow = NULL, margins = c(0.1, 24), 
              lwid = c(0.01, 4), lhei = c(2, 5), labCol = lC)
  }else{
    heatmap.2(keep.positions, trace = "none",
              breaks = 101, symbreaks = F, col = colors.hm, key = F,
              dendrogram = "column", Rowv = scar.hc.ro, Colv = celltree, 
              labRow = NULL, margins = c(0.1, 24), cexRow = scar.plot.size,
              lwid = c(0.01, 4), lhei = c(2, 5), labCol = lC)
  }
}
 
resample.links <- function(link.dist, fscars.r, n.sample){
  # Resample all links by redistributing the scars over the cells, calculate
  # the link ratio between properties, and add those to a distribution of
  # resampled link ratios. After all the resampling, calculate in which
  # percentile the original observation falls, and calculate the fold change
  # between original link ratio and mean resampled link ratios.
  property <- colnames(fscars.r)[3]
  
  max.scars.per.cell <- max(table(fscars.r$Cell))
  strata <- unlist(table(fscars.r[, 3]))
  
  scars <- fscars.r$Scar
  
  cells <- unique(fscars.r[, c(1, 3)])
  cells.to.sample <- data.frame(Cell = rep(cells$Cell, max.scars.per.cell))
  cells.to.sample <- merge(cells.to.sample, cells)
  sample.step.size <- ceiling(n.sample/10)
  # print("Starting sampling")
  
  for(s in 1:n.sample){
    #     print(paste("Sample", s))
    #     if((s %% sample.step.size) == 0){
    #       print(paste("Sampling step", s/sample.step.size))
    #     }
    #     # Shuffle scars
    sh.scars <- scars[sample.int(length(scars))]
    # Sample cells and multiplicities to assign scars to
    sampled.scars <- stratified(cells.to.sample, property, strata)
    # Paste scars to cells
    sampled.scars$Scar <- sh.scars
    # Up to this point, everything works hunky-dory, but it is possible for
    # cells to get the same scar twice. With 110 cells and 166 scars, some with
    # high multiplicities, this happens roughly all the time.
    # Switch the duplicate scars with scars in other cells and check again.
    while(sum(duplicated(sampled.scars)) > 0){
      # Take the first duplicated scar
      duplicated.row <- which(duplicated(sampled.scars))[1]
      # Select a random generated scar/cell from the same stratum (I don't
      # think this is necessary, and for strata with one loop this can easily
      # lead to an infinite loop)
      #       stratum <- 
      #         sampled.scars[sampled.scars[,2] == sampled.scars[duplicated.row, 2], ]
      #       replacement <- stratum[sample.int(nrow(stratum), 1), ]
      #       replacement.row <- 
      #         which(rownames(sampled.scars) == rownames(replacement))
      replacement.row <- sample.int(nrow(sampled.scars), 1)
      replacement <- sampled.scars[replacement.row, ]
      
      #       sampled.scars[duplicated.row, ]
      #       sampled.scars[replacement.row, ]
      #       # Switch the scars in both cells.
      duplicate.scar <- sampled.scars$Scar[duplicated.row]
      sampled.scars$Scar[duplicated.row] <- replacement$Scar  
      sampled.scars$Scar[replacement.row] <- duplicate.scar
      
      # sum(duplicated(sampled.scars))
      
    }
    
    sampled.scars <- sampled.scars[, c(1, 3, 2)]
    sampled.links <- calculate.links(sampled.scars)
    colnames(sampled.links)[3] <- paste("S", s, sep = ".")
    link.dist <- merge(link.dist, sampled.links)
  }
  
  link.dist$Depleted <- NA
  link.dist$Enriched <- NA
  link.dist$FC <- NA
  
  for(i in 1:nrow(link.dist)){
    x <- link.dist[i, 3:(n.sample + 3)]
    
    s <- as.numeric(x[-1])
    comp <- as.numeric(x[1])
    
    # Test for smaller than expected: how big a fraction of cases has a value that 
    # is smaller or equal to the observed value?
    link.dist$Depleted[i] <- sum(s <= comp)/length(s)
    
    # Test for larger than expected: how big a fraction of cases has a value that
    # is larger or equal to the observed value?
    link.dist$Enriched[i] <- sum(s >= comp)/length(s)
    
    # Calculate fold change between observation and randomized.
    link.dist$FC[i] <- comp/mean(s)
  }
  
  return(link.dist[, c(1:3, n.sample + 4:6)])
}

query.simulation <- function(simulation.results, vary, constants, plot = T){
  query.raw.output <- simulation.results
  plot.title <- "|"
  for(l in 1:length(constants)){
    selecter <- query.raw.output[, eval(names(constants)[l])] == constants[l]
    query.raw.output <- query.raw.output[selecter, ]
    plot.title <- 
      paste(plot.title, " ", names(constants)[l], ": ", constants[l], " |", sep = "")
  }
  query.raw.output <- query.raw.output[complete.cases(query.raw.output), ]
  query.aggregate.output <-
    aggregate(query.raw.output$true.pos.rate,
              by = list(p.A = query.raw.output$p.A,
                        expansions = query.raw.output$expansions,
                        Start.cells = query.raw.output$Start.cells),
              mean)
  colnames(query.aggregate.output)[4] <- "true.pos.rate"
  
  query.output.list <- list(raw = query.raw.output,
                            aggregate = query.aggregate.output)
  
  if(plot){
    title.elements <- paste(constants)
    print(
      ggplot() +
        geom_point(data = query.raw.output, size = 1,
                   aes_string(x = vary, y = "true.pos.rate")) +
        geom_point(data = query.aggregate.output, size = 2, color = "red",
                   aes_string(x = vary, y = "true.pos.rate")) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(title = plot.title,
             y = "True positive rate")
    )
  }
}

scar.tree.step <- function(scar.graph, scar.tree.in, 
                           scar.pref.order = character()){
  # Loop over connected components, remove scars in the preferred order or, if
  # none are present, determine the node(s) with the highest degree
  # and remove, write the name of the node to the dataframes, recursively call
  # function on remaining graph if there is one. Return a list of
  # edges and nodes created and the parent node, i.e. the current level of
  # recursion.
  if(length(scar.pref.order) == 0){
    scar.deletion.order <- data.frame(Scar = character(),
                                      Order = integer())
  }else{
    scar.deletion.order <- data.frame(Scar = scar.pref.order,
                                      Order = 1:length(scar.pref.order))
  }  
  decomposed.graph <- decompose(scar.graph)
  scar.tree.out <- scar.tree.in
  parent.node <- scar.tree.in$Parent.node
  
  # g <- 1 # For debugging
  for(g in 1:length(decomposed.graph)){
    scar.sub.graph <- decomposed.graph[[g]]
    scar.sub.graph.degrees <- 
      data.frame(Degree = degree(scar.sub.graph, mode = "all", loops = F))
    scar.sub.graph.degrees$Scar <- rownames(scar.sub.graph.degrees)
    scar.sub.graph.degrees <- 
      merge(scar.sub.graph.degrees, scar.deletion.order,
            all.x = T)
    scar.sub.graph.degrees$Order[is.na(scar.sub.graph.degrees$Order)] <-
      nrow(scar.deletion.order) + 1
    
    # If any of the scars present are on the scars preferred order list, delete
    # the one highest on that list. Otherwise, delete the the scar with the
    # highest degree.
    if(min(scar.sub.graph.degrees$Order) <= nrow(scar.deletion.order)){
      scar.to.delete <- 
        scar.sub.graph.degrees$Scar[scar.sub.graph.degrees$Order == 
                                      min(scar.sub.graph.degrees$Order)]
    }else{
      scar.to.delete <- 
        scar.sub.graph.degrees$Scar[scar.sub.graph.degrees$Degree == 
                                      max(scar.sub.graph.degrees$Degree)]
    }
    # print(cat("Deleting scar", scar.to.delete))
    # if(nrow(scar.tree.out$Nodes) == 0){ 
    #   current.node <- parent.node + 1
    # }else{
    #   current.node <- max(scar.tree.out$Nodes$Node.number) + 1
    current.node <- paste(parent.node, g, sep = "_")
    # }
    
    scar.tree.out$Nodes <-
      rbind(scar.tree.out$Nodes,
            data.frame(Node = current.node,
                       Scar = scar.to.delete))
    # print(scar.tree.out$Node)
    scar.tree.out$Edges <- 
      rbind(scar.tree.out$Edges,
            data.frame(From = parent.node,
                       To = current.node))
    # print(scar.tree.out$Edges)
    
    scar.graph.down <- delete.vertices(scar.sub.graph, scar.to.delete)
    
    continue.downward <- length(V(scar.graph.down)) > 0
    
    # parent.node <- current.node # For debugging
    # scar.graph <- scar.graph.down # For debugging
    if(continue.downward){
      scar.tree.out$Parent.node <- current.node
      scar.tree.out <- scar.tree.step(delete.vertices(scar.sub.graph, scar.to.delete),
                                      scar.tree.out, scar.pref.order)
    }
  }
  
  return(scar.tree.out)
}

scar.tree.step.old <- function(scar.graph, scar.tree.in){
  # Loop over connected components, determine the node(s) with the highest degree
  # and remove, write the name of the node to the dataframes, recursively call
  # function on remaining graph if there is one. Return a list of
  # edges and nodes created and the parent node, i.e. the current level of
  # recursion.
  decomposed.graph <- decompose(scar.graph)
  scar.tree.out <- scar.tree.in
  parent.node <- scar.tree.in$Parent.node
  
  # g <- 1 # For debugging
  for(g in 1:length(decomposed.graph)){
    scar.sub.graph <- decomposed.graph[[g]]
    scar.sub.graph.degrees <- 
      data.frame(Degree = degree(scar.sub.graph, mode = "all", loops = F))
    scar.sub.graph.degrees$Scar <- rownames(scar.sub.graph.degrees)
    
    scar.to.delete <- 
      scar.sub.graph.degrees$Scar[scar.sub.graph.degrees$Degree == 
                                    max(scar.sub.graph.degrees$Degree)]
    # print(cat("Deleting scar", scar.to.delete))
    # if(nrow(scar.tree.out$Nodes) == 0){ 
    #   current.node <- parent.node + 1
    # }else{
    #   current.node <- max(scar.tree.out$Nodes$Node.number) + 1
    current.node <- paste(parent.node, g, sep = "_")
    # }
    
    scar.tree.out$Nodes <-
      rbind(scar.tree.out$Nodes,
            data.frame(Node = current.node,
                       Scar = scar.to.delete))
    # print(scar.tree.out$Node)
    scar.tree.out$Edges <- 
      rbind(scar.tree.out$Edges,
            data.frame(From = parent.node,
                       To = current.node))
    # print(scar.tree.out$Edges)
    
    scar.graph.down <- delete.vertices(scar.sub.graph, scar.to.delete)
    
    continue.downward <- length(V(scar.graph.down)) > 0
    
    # parent.node <- current.node # For debugging
    # scar.graph <- scar.graph.down # For debugging
    if(continue.downward){
      scar.tree.out$Parent.node <- current.node
      scar.tree.out <- scar.tree.step(delete.vertices(scar.sub.graph, scar.to.delete),
                                      scar.tree.out)
    }
  }
  
  return(scar.tree.out)
}


simulate.p.right <- function(p.A, expansion, n.c){
  # Calculate statistical chance of calling a link correct by simulation.
  
  # Simulate the number of scars created.
  scar.numbers <- rbinom(1, n.c, p.A)
  
  total.connections <- choose(n.c * 2^expansion, 2)
  
  # The total number of connections between a pool of cells that all have scar
  # A. For a pool of size [scar.A.created] * 2^[expansion], this is
  # choose([scar.A.created] * 2^[expansion], 2)
  total.scar.connections <- choose(scar.numbers * 2^expansion, 2)
  # The connections between cells with scar A that are related through 
  # expansion
  real.connections <- scar.numbers * choose(2^expansion, 2)
  # The connections between cells with scar A that are unrelated through
  # expansion. Size of every clone (there are [scar.A.created] clones) 
  # is 2^[expansion], number of connections between two clones then is
  # 2^(2 * [expansion]). There are choose([scar.A.created], 2) combinations
  # between clones to consider.
  false.connections <- 2^(2 * expansion) * choose(scar.numbers, 2)
  
  t.p.r <- real.connections/total.scar.connections
  f.p.r <- false.connections/total.scar.connections

  return(c(t.p.r, f.p.r))
}

simulate.p.right.2 <- function(p.A, expansion, n.c, trials){
  # Calculate chance of calling a link correct by simulation.
  scar.numbers.exp <- data.frame(Trial = 1:trials,
                                 scar.A.created = rbinom(trials, n.c, p.A))
  # The total number of connections possible between a cell pool that started
  # with [n.c] cells that each expanded [expansion] times.
  scar.numbers.exp$Total.connections <- choose(n.c * 2^expansion, 2) 
  # The total number of connections between a pool of cells that all have scar
  # A. For a pool of size [scar.A.created] * 2^[expansion], this is
  # choose([scar.A.created] * 2^[expansion], 2)
  scar.numbers.exp$Total.A.connections <-
    choose(scar.numbers.exp$scar.A.created * 2^expansion, 2)
  # The connections between cells with scar A that are related through 
  # expansion
  scar.numbers.exp$Real.connections <- 
    scar.numbers.exp$scar.A.created * choose(2^expansion, 2)
  # The connections between cells with scar A that are unrelated through
  # expansion. Size of every clone (there are [scar.A.created] clones) 
  # is 2^[expansion], number of connections between two clones then is
  # 2^(2 * [expansion]). There are choose([scar.A.created], 2) combinations
  # between clones to consider.
  scar.numbers.exp$False.connections <-
    2^(2 * expansion) * choose(scar.numbers.exp$scar.A.created, 2)
  scar.numbers.exp$p.right <- 
    scar.numbers.exp$Real.connections/scar.numbers.exp$Total.A.connections
  
  return(scar.numbers.exp)
}
