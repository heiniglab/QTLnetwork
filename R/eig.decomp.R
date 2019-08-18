#' -----------------------------------------------------------------------------
#' Get eigenvalue decomposition
#'
#' @param M
#' @param n.eigs
#' @param sym
#'
#' @return
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
eig.decomp <- function(M, n.eigs, sym) {

  n = nrow(M)
  f <- function(x, extra=NULL) as.matrix(M %*% x)
  wh <- if (sym) 'LA' else 'LM'

  # constraints: n >= ncv > nev
  ar <- arpack(f, sym = sym, options = list(
               which = wh, n = n, ncv = min(n, 4*n.eigs), nev = n.eigs + 1))
  if (!sym) {
    ar$vectors <- Re(ar$vectors)
    ar$values  <- Re(ar$values)
  }

  # check for negative values and flip signs
  neg = which(ar$values < 0)
  for (n in neg) {
    ar$values[n] = -ar$values[n]
    ar$vectors[,n] = -ar$vectors[,n]
  }

  return(ar)
}
