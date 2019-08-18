#' -----------------------------------------------------------------------------
#' Get nodes by type
#'
#' Convenience function to retrieve nodes of a certain type from a graph
#' object using nodeData attributes
#'
#' @param graph graphNEL (locus graph) object
#' @param type character indicating the node type
#'
#'
#' @return Character vector containing nodes corresponding to the given type
#'
#' @import graph
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
get.nodes.by.type <- function(graph, type) {
  require(graph)
  n = nodes(graph)
  return(n[unlist(nodeData(graph, n, type))])
}
