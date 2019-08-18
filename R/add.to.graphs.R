#' -----------------------------------------------------------------------------
#' Augment grahs
#'
#' Add TF to CpG edges and trans gene to SNP edges to graphs. This function
#' assumes that you have an object tfbs.ann in the environment. It is not
#' passed as argument (bad style) becuase it is too big (passing by value).
#'
#' @param graphs list of graph objects to add nodes and edges to
#' @param sentinel character id of the sentinel SNP
#' @param trans.genes character vector of genes in the trans locus
#' @param trans.cpgs character vector of CpGs with trans meQTLs
#' @param tfbs.ann logical indicator matrix with CpGs in the rows and
#'        transcription factors as columns
#'
#' @return list of graph objects with the nodes and edges added
#'
#' @import reshape
#' @import graph
#'
#' @export
#'
#' @author Matthias Heinig <matthias.heinig@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
add.to.graphs <- function(graphs,
                          sentinel,
                          trans.genes,
                          trans.cpgs,
                          tfbs.ann) {
  require(reshape)
  require(graph)

  tf.edges = melt(tfbs.ann[trans.cpgs,])
  colnames(tf.edges) = c("cpg", "condition", "adjacent")
  tf.edges = tf.edges[tf.edges[,"adjacent"],]
  for (i in 1:2) {
    tf.edges[,i] = as.character(tf.edges[,i])
  }
  tf = as.character(sapply(strsplit(tf.edges[,"condition"], ".", fixed=T),
                           "[", 1))
  tf.edges = data.frame(tf.edges, tf, stringsAsFactors=F)
  tf.edges = tf.edges[!duplicated(paste(tf.edges[,"tf"],
                                        tf.edges[,"cpg"])),]

  ## put together the graphs
  out = list()

  for (graph.idx in 1:length(graphs)) {
    locus.graph = graphs[[graph.idx]]

    ## filter for TFs that are in the graph already
    use.tf.edges = tf.edges[tf.edges[,"tf"] %in% nodes(locus.graph),]

    new.nodes = unique(c(use.tf.edges[,"tf"],
                         trans.cpgs,
                         sentinel,
                         trans.genes))
    new.nodes = setdiff(new.nodes, nodes(locus.graph))

    locus.graph = addNode(new.nodes, locus.graph)

    ## also add some meta data (the type of nodes)
    nodeDataDefaults(locus.graph, "tf") = FALSE
    nodeData(locus.graph, unique(use.tf.edges[,"tf"]), "tf") = TRUE

    nodeDataDefaults(locus.graph, "cpg") = FALSE
    nodeData(locus.graph, trans.cpgs, "cpg") = TRUE

    nodeDataDefaults(locus.graph, "snp") = FALSE
    nodeData(locus.graph, sentinel, "snp") = TRUE

    nodeDataDefaults(locus.graph, "trans.gene") = FALSE
    nodeData(locus.graph, trans.genes, "trans.gene") = TRUE

    ## add edges for the tfbs
    locus.graph = addEdge(use.tf.edges[,"tf"],
                          use.tf.edges[,"cpg"],
                          locus.graph)

    ## add edges for the connection of the locus and its genes
    locus.graph = addEdge(rep(sentinel, length(trans.genes)),
                          trans.genes,
                          locus.graph)

    out[[graph.idx]] = locus.graph
  }

  return(out)
}
