#' -----------------------------------------------------------------------------
#' Convert graph to a sparse matrix
#'
#' @param graph The graph to be converted to sparse matrix representation
#'
#' @return The graph represented as a sparse matrix
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
graph2sparseMatrix <- function(graph) {
  em = edgeMatrix(graph)
  Asparse = sparseMatrix(em[1,], em[2,], x=1, dims=rep(numNodes(graph), 2), dimnames=list(nodes(graph), nodes(graph)))

  ## make symmetric
  Asparse = Asparse + t(Asparse)
  Asparse[Asparse > 1] = 1

  return(Asparse)
}
