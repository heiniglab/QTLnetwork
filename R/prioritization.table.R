#' -----------------------------------------------------------------------------
#' Prioritization summary table
#'
#' Extract summary information from the random walk analysis
#'
#' @param prioritize list of random walk analysis results generated with the
#'                  network_batch.R script
#' @param full boolean indicating whether all genes in the trans locus should
#'             be reported (TRUE) or only the best (FALSE)
#' @param locus.graphs a list of graphNEL objects with more details on the
#'                     graphs used for analysis (for full summary)
#' @param subset list of subsets of trans genes for each locus this is useful
#'               to recycle permutation runs on maximal intervals
#'
#' @return data.frame with the best candidate gene per locus
#'
#' @import graph
#'
#' @export
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
prioritization.table <- function(prioritize, full=FALSE,
                                 locus.graphs=NULL, subset=NULL) {
  require(graph)
  tab = NULL
  for (sentinel in names(prioritize)) {
    if (length(prioritize[[sentinel]]$best.trans.gene) == 0) {
      cat("problem with SNP:", sentinel, "\n")
      next
    }

    if (full) {
      if (is.null(locus.graphs)) {
        stop("need to specify the locus.graphs list when full == TRUE")
      }

      ## get the full list of trans genes (also the ones not in string)
      all.tgenes = get.nodes.by.type(locus.graphs[[sentinel]], "trans.gene")
    } else {
      all.tgenes = with(prioritize[[sentinel]], trans.genes)
    }
    ## check if we are considering only a subset
    if (!is.null(subset)) {
      sid = strsplit(sentinel, ".", fixed=T)[[1]][1]
      ## if (sid == "rs7783715") {
      ##   browser()
      ## }
      all.tgenes = intersect(all.tgenes, subset[[sid]])
    }

    ## put together the details for the tested trans genes
    tested.tgenes = intersect(with(prioritize[[sentinel]], trans.genes),
      all.tgenes)
    if (length(tested.tgenes) == 0) {
      cat("No tested trans genes for SNP:", sentinel, "\n")
      next
    }
    best.tgene = with(prioritize[[sentinel]],
      tested.tgenes[which.max(prop[tested.tgenes,"from"])])

    p = with(prioritize[[sentinel]], prop[tested.tgenes,"from"])
    rescaled = p / sum(p)
    entropy = sum(-rescaled * log(rescaled))
    ## also compute the maximum possible
    u = rep(1/length(p), length(p))
    max.entropy = sum(-u * log(u))

    trans.gene.details = with(prioritize[[sentinel]],
      data.frame(trans.gene=tested.tgenes, prop=prop[tested.tgenes,"from"],
                 p, in.string=TRUE))

    ## add the genes that were not in string
    not.in.string = setdiff(all.tgenes, tested.tgenes)
    if (length(not.in.string) > 0) {
      trans.gene.details = rbind(trans.gene.details,
        data.frame(trans.gene=not.in.string, prop=NA, p=NA, in.string=FALSE))
    }

    smry = with(prioritize[[sentinel]],
      data.frame(sentinel, ncpg=length(cpgs), ntrans=length(tested.tgenes),
                 best.trans=best.tgene, best.tf=best.tf,
                 trans.entropy=entropy, max.entropy=max.entropy,
                 rel.entropy=entropy/max.entropy,
                 best.trans.p=max(rescaled), best.trans.prop=max(p),
                 trans.gene.details))

    ## just get the best trans gene
    if (!full) {
      smry = smry[smry[,"trans.gene"] == best.tgene,]
    }

    tab = rbind(tab, smry)
  }
  return(tab)
}
