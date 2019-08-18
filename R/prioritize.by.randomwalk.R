#' -----------------------------------------------------------------------------
#' Candidate prioritization by random walk analysis
#'
#' @param locus.graph graphNEL object as generated in network_preprocess.R
#' @param largest.cc graphNEL object with the largest connected component
#'        of the protein - protein interaction network
#' @param string.nodes character vector with the names of the nodes in the PPI
#'        network
#' @param tf.nodes character vector with the names of the transcription factors
#' @param n.eigs number of eigen vectors to use for the approximation
#' @param mode character "propagation" or "hitting.time"
#' @param return.g logical should the graph used for analysis be returned?
#'
#' @return list with sentinel, cpgs, tfs, trans.genes, prop,
#'   best.trans.gene, best.tf, g
#'
#' @import graph
#' @import igraph
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
prioritze.by.randomwalk <- function(locus.graph, largest.cc, string.nodes,
                                    tf.nodes, n.eigs=500, mode="propagation",
                                    return.g=FALSE) {

  require(graph)
  require(igraph)

  n = nodes(locus.graph)

  ## get the sentinel
  sentinel = n[unlist(nodeData(locus.graph, n, "snp"))]

  ## get all the trans locus genes
  trans.genes = n[unlist(nodeData(locus.graph, n, "trans.gene"))]
  trans.genes = intersect(trans.genes, graph::nodes(largest.cc))
  if (length(trans.genes) == 0) {
    warning("no trans genes found in locus graph")
    return(NULL)
  }

  ## also get the cpgs
  cpgs = n[unlist(nodeData(locus.graph, n, "cpg"))]
  ## cpgs = cpgs[grep("cg", cpgs)]
  cpgs = intersect(cpgs, graph::nodes(largest.cc))
  if (length(cpgs) == 0) {
    warning("no cpgs found in locus graph")
    return(NULL)
  }

  ## get the sub graph for the locus (removing snps and all other cpgs)
  ## KAP1 is not connected in the STRING network
  nodes.set = c(string.nodes, setdiff(tf.nodes, "KAP1"), trans.genes, cpgs)
  g = subGraph(intersect(nodes(largest.cc), nodes.set), largest.cc)

  ## check again that we have one connected component
  ig = graph_from_graphnel(g)
  cl = clusters(ig)
  keep = nodes(g)[cl$membership == which.max(cl$csize)]
  if (!all(nodes(g) %in% keep)) {
    missing = setdiff(nodes(g), keep)
    cat("Some nodes are not in the connected component! Please check input!\n")
    cat(missing, sep="\n")
    cat("\n")
    missing = cbind(absolute=sapply(list(trans.genes=trans.genes,
                      string=string.nodes, tfs=tf.nodes), function(x)
                      length(intersect(x, missing))),
      relative=sapply(list(trans.genes=trans.genes, string=string.nodes,
        tfs=tf.nodes), function(x)
        length(intersect(x, missing)) / length(x)))
    print(missing)
    g = subGraph(keep, largest.cc)
    cpgs = intersect(cpgs, nodes(g))
    trans.genes = intersect(trans.genes, nodes(g))
    tf.nodes = intersect(tf.nodes, nodes(g))
  }

  rg = NULL
  if (return.g) {
    rg = g
  }

  ## also get the tfs
  tfs = unique(unlist(adj(g, cpgs)))

  ## now we would like to get the propagation values from the cpgs to the
  ## genes in the trans locus and the tfs
  if (mode == "propagation") {
    prop = propagation(graph2sparseMatrix(g), n.eigs=n.eigs,
                       from=cpgs, to=c(trans.genes, tfs), sum="both")

    best.trans = trans.genes[which.max(prop[trans.genes,"from"])]

    return(list(sentinel=sentinel, cpgs=cpgs, tfs=tfs, trans.genes=trans.genes,
                prop=prop, best.trans.gene=best.trans,
                best.tf=tfs[which.max(prop[tfs,"from"])], g=rg))
  }

  if (opt$mode == "hitting.time") {
    ht = hitting.time(graph2sparseMatrix(g), n.eigs=500, from=cpgs,
                      to=c(trans.genes, tfs))
    return(list(sentinel=sentinel,
                cpgs=cpgs,
                tfs=tfs,
                trans.genes=trans.genes,
                ht=ht,
                best.trans.gene=trans.genes[which.max(colSums(ht[,trans.genes,drop=F]))],
                best.tf=tfs[which.max(colSums(ht[,tfs,drop=F]))], g=rg))
  }

}
