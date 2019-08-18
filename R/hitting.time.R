#' -----------------------------------------------------------------------------
#' Alternative approach to prioritize genes
#'
#' We can compute expected hitting times for the random walks instead of using
#' the propagation approach.
#'
#' @details Using \code{from} and \code{to} nodes avoids to compute the full
#'   propagation matrix which is too large to fit in memory. When \code{from} and
#'   \code{to} are NULL the full matrix will be computed. For \code{sum}, "none"
#'   returns a "from" by "to" matrix; "from" returns a vector of length
#'   nrow(Asparse); "to" returns a vector of length nrow(Asparse) and "both"
#'   returns a nrow(Asparse) by 2 matrix with the \code{from} and \code{to}
#'   vectors. The i-th component of the "from sum" represents the average
#'   hitting time going from node i to any of the "to" nodes. The i-th component
#'   of the "to sum" represents the average hitting time going from any "from"
#'   node to node i
#'
#' @param n.eigs The number of eigenvectors used for the approximation. If this
#'   is smaller than 2, the pseudo inverse will be computed
#' @param from Source nodes for which we need the transition probabilties.
#' @param to Target nodes for which we need the transition probabilties.#'
#' @param sum Can be "none", "from", "to" or "both" and will sum up the
#'   propagation matrix over one or the other or both lists.
#'
#' @return Table of hitting times
#'
#' @export
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
hitting.time <- function(Asparse, n.eigs=20, from=NULL, to=NULL, sum="none") {

  ## we need the number of edges and the number of nodes and the degree of nodes
  numEdges = sum(Asparse) / 2 + sum(diag(Asparse))
  numNodes = nrow(Asparse)
  deg = rowSums(Asparse)

  ## symmetric transition matrix
  D.root = Diagonal(nrow(Asparse), 1 / sqrt(rowSums(Asparse)))
  Psym = D.root %*% Asparse %*% D.root

  if (is.null(from)) {
    from = 1:numNodes
  }
  if (is.null(to)) {
    to = 1:numNodes
  }
  if (is.logical(from)) {
    from = which(from)
  }
  if (is.logical(to)) {
    to = which(to)
  }
  dnames = list(from, to)
  ## if we have character we need to match
  if (is.character(from)) {
    from = match(from, rownames(Asparse))
  }
  if (is.character(to)) {
    to = match(to, colnames(Asparse))
  }

  ## setup different ways of summarizing the hitting.time matrix
  if (sum == "from") {
    transform = rowMeans
    hitting.time = double(0, length(to))
    names(hitting.time) = dnames[[2]]
  } else if (sum == "to") {
    transform = colMeans
    hitting.time = double(0, length(from))
    names(hitting.time) = dnames[[1]]
  } else if (sum == "both") {
    hitting.time = matrix(0, nrow=numNodes, ncol=2,
                          dimnames=list(rownames(Asparse), c("from", "to")))
  } else if (sum == "none") {
    transform = function(x) {return(x)}
    hitting.time = matrix(0, nrow=length(from), ncol=length(to), dimnames=dnames)
  }

  ## approximate using the eigen vectors

  if (n.eigs < nrow(Psym)) {
    ## just use a subset of eigen vectors
    eig.M <- eig.decomp(Psym, n.eigs, TRUE)

  } else {
    ## use all eigen vectors (only for small matrices)
    eig.M <- eigen(Psym)
  }

  ## Hitting time computed according to Theorem 3.1 from
  ## Lov??sz, L. (1993). Random walks on graphs. Combinatorics.

  compute.hitting.time <- function(from.nodes, to.nodes) {
    sapply(to.nodes, function(tt)
           sapply(from.nodes, function (ff) {
             ht = 2 * numEdges * sum(sapply(2:n.eigs, function(i) {
               v = eig.M$vectors[,i]
               return(1 / (1 - eig.M$values[i]) *
                      (v[tt]^2 / deg[tt] - v[ff] *
                       v[tt] / sqrt(deg[ff] * deg[tt])))
             }))
           }))
  }

  if (sum %in% c("none", "from", "to")) {
    for (i in 2:n.eigs) {
      increment = compute.hitting.time(from, to)
      hitting.time = hitting.time + transform(increment)
    }

  } else {
    ## increment the "from" sum vector
    for (tt in to) {
      hitting.time[,1] = hitting.time[,1] + compute.hitting.time(1:numNodes, tt)
    }
    hitting.time[,1] = hitting.time[,1] / length(to)

    ## increment the "to" sum vector
    for (ff in from) {
      hitting.time[,2] = hitting.time[,2] + compute.hitting.time(ff, 1:numNodes)
    }
    hitting.time[,2] = hitting.time[,2] / length(from)
  }

  return(hitting.time)
}
