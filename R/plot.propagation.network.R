#' -----------------------------------------------------------------------------
#' Convenience function to plot networks
#'
#'
#' @import graph
#' @import Rgraphviz
#' @import RColorBrewer
#'
#' @export
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
plot.propagation.network <- function(g, prop, cpgs, best.trans, trans.genes,
                                     sentinel, prefix, additional.genes=NULL,
                                     additional.bordercol=NULL) {
  require(graph)
  require(Rgraphviz)
  require(RColorBrewer)

  if (!file.exists(dirname(prefix))) {
    dir.create(dirname(prefix), recursive=T)
  }

  ## now we get the node weights by adding up the from and to weights
  node.weight = rowSums(prop)

  ## finally we would like to find the shortest path with maximal weight
  ## we have an algorithm that finds minimum node weight paths so we need
  ## to turn the weighting around
  ## in addition weights need to be non-negative
  node.weight = max(node.weight) - node.weight + 1

  sp = min.node.weight.path(g, node.weight, from=cpgs, to=best.trans)

  plot.nodes = unique(c(setdiff(unlist(lapply(sp, "[", "path_detail")), NA),
                        sentinel, trans.genes, additional.genes))

  plot.edges = NULL
  for (i in sp) {
    if (length(i$path_detail) > 0) {
      for (j in 1:(length(i$path_detail) - 1)) {
        plot.edges = rbind(plot.edges,
          data.frame(from=i$path_detail[j], to=i$path_detail[j + 1],
                     stringsAsFactors=F))
      }
    }
  }

  ## plot.graph = graphNEL(nodes=unique(c(plot.nodes, trans.genes)))
  ## plot.graph = addEdge(sentinel, trans.genes, plot.graph)
  ## plot.graph = addEdge(plot.edges[,1], plot.edges[,2], plot.graph)
  plot.graph = g
  plot.graph = graph::addNode(sentinel, plot.graph)
  plot.graph = graph::addEdge(sentinel, trans.genes, plot.graph)
  plot.graph = graph::subGraph(plot.nodes, plot.graph)

  ## instead of selecting just one shortest path, we could also use a modified
  ## Dijkstra algorithm that finds all shortest paths

  attrs <- list(node=list(fixedsize=TRUE, fontsize=10, style="filled",
                          fontname="helvetica"),
                graph=list(overlap="true",
                           root=sentinel,
                           outputorder="edgesfirst"))

  shape = rep("ellipse", numNodes(plot.graph))
  names(shape) = nodes(plot.graph)
  shape[grep("cg", nodes(plot.graph))] = "box"
  shape[sentinel] = "box"

  width = rep(0.8, numNodes(plot.graph))
  names(width) = nodes(plot.graph)
  width[grep("cg", nodes(plot.graph))] = 0.2

  height = rep(0.2, numNodes(plot.graph))
  names(height) = nodes(plot.graph)
  height[grep("cg", nodes(plot.graph))] = 0.2

  label = nodes(plot.graph)
  names(label) = nodes(plot.graph)
  label[grep("cg", nodes(plot.graph))] = ""

  ## fill color represents the random walk score
  pal = brewer.pal(9, "Blues")
  score = prop[match(nodes(plot.graph), rownames(prop)), "from"]
  names(score) = nodes(plot.graph)
  breaks = seq(min(score, na.rm=T), max(score, na.rm=T),
               length.out=length(pal)+1)
  col = pal[findInterval(score, breaks, rightmost.closed=T)]
  names(col) = nodes(plot.graph)
  col[sentinel] = "#ff0000"
  col[grep("cg", nodes(plot.graph))] = "#00ff00"

  penwidth = rep(1, numNodes(plot.graph))
  names(penwidth) = nodes(plot.graph)
  penwidth[best.trans] = 2

  bordercol = rep("black", numNodes(plot.graph))
  names(bordercol) = nodes(plot.graph)
  bordercol[best.trans] = "red"

  if (!(is.null(additional.genes) || is.null(additional.bordercol))) {
    bordercol[additional.genes] = additional.bordercol
    penwidth[additional.genes] = 2
  }

  nAttrs = list(shape=shape, label=label, width=width, height=height,
                fillcolor=col, penwidth=penwidth, color=bordercol)

  ecol = rep("black", numEdges(plot.graph))
  names(ecol) = edgeNames(plot.graph)
  ecol[grep(sentinel, names(ecol))] = "red"

  eAttrs = list(color=ecol)


  pdf.file = paste(prefix, ".pdf", sep="")
  dot.file = paste(prefix, ".dot", sep="")
  rdata.file = paste(prefix, ".RData", sep="")

  ## direct plotting plots edges over nodes
  pdf(file=pdf.file)
  plot(plot.graph, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  dev.off()

  ## also use the graphviz command line for plotting
  toDot(plot.graph, dot.file, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  system(paste("twopi -Tpdf -O", dot.file))

  ## plot the color legend
  pdf(file=paste(prefix, "-color-legend.pdf", sep=""))
  plot(c(0, 1), c(0, 9))
  rect(rep(0, 9), 0:8, rep(1, 9), 1:9, col=pal)
  dev.off()

  ## also save the data for later
  save(list=c("nAttrs", "eAttrs", "attrs", "plot.graph"), file=rdata.file)

  ## browser()
}
