#' -----------------------------------------------------------------------------
#' Approximate the infinite sum
#'
#' In the paper and notes this is matrix is called M#'
#'
#' @details Using \code{from} and \code{to} nodes avoids to compute the full
#'   propagation matrix which is too large to fit in memory. When \code{from} and
#'   \code{to} are NULL the full matrix will be computed. For \code{sum}, "none"
#'   returns a "from" by "to" matrix; "from" returns a vector of length
#'   nrow(Asparse); "to" returns a vector of length nrow(Asparse) and "both"
#'   returns a nrow(Asparse) by 2 matrix with the \code{from} and \code{to}
#'   vectors.
#'
#' @param n.eigs The number of eigenvectors used for the approximation. If this
#'   is smaller than 2, the pseudo inverse will be computed
#' @param from Source nodes for which we need the transition probabilties.
#' @param to Target nodes for which we need the transition probabilties.#'
#' @param sum Can be "none", "from", "to" or "both" and will sum up the
#'   propagation matrix over one or the other or both lists.
#'
#' @import igraph
#' @import Matrix
#' 
#' @export
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
propagation <- function(Asparse, n.eigs=20, from=NULL, to=NULL, sum="none") {

  require(igraph)
  require(Matrix)
  ## transition matrix
  ## transition = Asparse / rowSums(Asparse)

  ## symmetric transition matrix
  D.root = Diagonal(nrow(Asparse), 1 / sqrt(Matrix::rowSums(Asparse)))
  Psym = D.root %*% Asparse %*% D.root

  if (is.null(from)) {
    from = 1:nrow(Asparse)
  }
  if (is.null(to)) {
    to = 1:nrow(Asparse)
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

  ## setup different ways of summarizing the propagation matrix
  if (sum == "from") {
    transform = rowSums
    propagation = double(0, length(to))
    names(propagation) = dnames[[2]]
  } else if (sum == "to") {
    transform = colSums
    propagation = double(0, length(from))
    names(propagation) = dnames[[1]]
  } else if (sum == "both") {
    propagation = matrix(0, nrow=nrow(Asparse), ncol=2,
                         dimnames=list(rownames(Asparse), c("from", "to")))
  } else if (sum == "none") {
    transform = function(x) {return(x)}
    propagation = matrix(0, nrow=length(from), ncol=length(to), dimnames=dnames)
  }

  if (n.eigs > 0) {
    ## approximate using the eigen vectors

    if (n.eigs < nrow(Psym)) {
      ## just use a subset of eigen vectors
      eig.M <- eig.decomp(Psym, n.eigs, TRUE)

    } else {
      ## use all eigen vectors (only for small matrices)
      eig.M <- eigen(Psym)
    }

    ## for both computations we discard the first eigenvector as it
    ## represents the stationary distribution
    if (sum %in% c("none", "from", "to")) {
      for (i in 2:n.eigs) {
        increment = (eig.M$values[i] / (1 - eig.M$values[i])) *
          eig.M$vectors[from,i] %*% t(eig.M$vectors[to,i])
        propagation = propagation + transform(increment)
      }

    } else {
      ## when sum is "both" we iterate over the "from" and "to" nodes and add
      ## the from and to vectors incrementally to save memory
      for (i in 2:n.eigs) {
        weight = (eig.M$values[i] / (1 - eig.M$values[i]))
        ## increment the "from" sum vector
        for (ff in from) {
          propagation[,1] = propagation[,1] + weight *
            (eig.M$vectors[ff,i] * eig.M$vectors[,i])
        }
        ## increment the "to" sum vector
        for (tt in to) {
          propagation[,2] = propagation[,2] + weight *
            (eig.M$vectors[tt,i] * eig.M$vectors[,i])
        }
      }
    }
  } else {
    ## compute using the pseudo inverse

    ## the first eigenvector corresponds to stationary distribution defined
    ## by the degrees of the nodes
    ## Attention: this is actually the first eigenvector of the asymmetric
    ## Psym matrix
    d = rowSums(Asparse)
    v = sum(d)
    phi0 = d / v

    ## there is an equivalence of the eigenvectors of the symmetric and
    ## asymetric matrix:
    ## EV(sym) = D^-0.5 EV(asym)
    phi0 =  phi0 / sqrt(d)

    ## in the dtp package there is a length normalization step that matches
    ## the vectors exactly (normalizing to unit length vectors)
    phi0 = phi0 / sqrt(sum(phi0^2))

    n <- nrow(Psym)
    inv <- solve(Diagonal(n) - Psym + phi0 %*% t(phi0))
    propagation = inv - Diagonal(n)
    propagation = propagation[from, to]
    if (sum %in% c("none", "from", "to")) {
      propagation = transform(propagation)
    } else {
      stop("sum = 'both' not implemented yet for pseudo inverse")
    }
  }
  return(propagation)
}
