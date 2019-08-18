#' -----------------------------------------------------------------------------
#' Find the shortest path with minimal node weight.
#'
#' @param g The graphNEL object
#' @param weights The weights calculated for the nodes
#' @param from The start nodes from which to start the min-node-weight path
#' @param to The target nodes to which to calcaulate the min-node-weight path
#'
#' @import graph
#' @import RBGL
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
min.node.weight.path <- function(g, weights, from, to) {
  ## the problem can be transformed into a directed graph problem where all
  ## incoming edges are assigned the node weight

  require(graph)
  require(RBGL)

  h = graphNEL(nodes(g), edgemode="directed")
  em = edgeMatrix(g)
  h = addEdge(nodes(g)[em[1,]], nodes(g)[em[2,]], h, weights[em[2,]])
  h = addEdge(nodes(g)[em[2,]], nodes(g)[em[1,]], h, weights[em[1,]])

  return(sp.between(h, from, to))
}
