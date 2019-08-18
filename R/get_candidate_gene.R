#' -----------------------------------------------------------------------------
#' Get the candidate gene of a locus using the random walk approach.
#'
#' Applys the random walk approach to get the most likely SNP gene
#' from the list of provided SNP genes
#'
#' @param snp_id Character id of the SNP in the locus
#' @param snp_genes Character vector of the genes encoded near the SNP
#' @param trans_cpgs Character vector of the trans CpGs associated
#'   to the SNP
#' @param ppi_graph graphNEL object containing the PPIs (e.g. from STRING)
#' @param tfbs_annotation Logical indicator matrix with trans entities 
#'   in the rows and transcription factors as columns (bincinding site
#'   information of transcription factors at trans entities)
#' @param gene_only Logical whether to return only the identified candidate
#'   gene or also the rest of the random walk results (default: FALSE)
#'   
#' @return The candidate gene as determined by the random walk analysis
#'
#' @export
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#'#'@examples
#'\dontrun{
#' data(ppi_string)
#' data(test_data)
#' data(tfbs_annotation)
#' 
#' candidate_gene <- get_candidate_gene(snp_id, 
#'                                      snp_genes, 
#'                                      trans_cpgs, 
#'                                      ppi_string,
#'                                      tfbs_annotation,
#'                                      gene_only = TRUE)
#' print(candidate_gene)
#'}
#' -----------------------------------------------------------------------------
get_candidate_gene <- function(snp_id, 
                               snp_genes, 
                               trans_cpgs, 
                               ppi_graph,
                               tfbs_annotation,
                               gene_only = FALSE) {
  require(igraph)
  require(graph)
    
  # add trans genes and SNP genes to graph if needed. Use TFBS annotation
  # to connect trans 
  locus_graph <- add.to.graphs(list(ppi_graph),
                        snp_id, 
                        snp_genes, 
                        trans_cpgs, 
                        tfbs_annotation)[[1]]
  
  print("Creating graphs.")
  
  # this was used to speed up processing for a large number of loci
  # in the permutation analysis. We don't really need this individual
  # graph for a single locus, but it is needed in prioritize.by.randomwalk
  # at the moment
  joined_graph_ppi <- locus_graph
  ig = graph_from_graphnel(joined_graph_ppi)
  cl = clusters(ig)
  keep = nodes(joined_graph_ppi)[cl$membership == which.max(cl$csize)]
  joined_graph_ppi = subGraph(keep, joined_graph_ppi)
  
  ## gather relevant nodes sets for the random walk
  snp_genes_in_ppi <- snp_genes[snp_genes %in% nodes(ppi_graph)]
  cpgs_with_tfbs <- trans_cpgs[graph::degree(locus_graph, trans_cpgs) > 0]
  # we expect colnames to be of the format 'TF.cell-type'
  tf_nodes = setdiff(intersect(nodes(joined_graph_ppi), 
                       sapply(strsplit(colnames(tfbs_annotation), ".", fixed=T), 
                              "[", 1)), "KAP1")
  nodeset = unique(c(nodes(ppi_graph), snp_genes_in_ppi, 
                     cpgs_with_tfbs, snp_id))
  joined_graph_ppi = subGraph(intersect(nodes(joined_graph_ppi), nodeset), joined_graph_ppi)
  
  print("Applying random walk.")
  # apply the random walk priorization, this is prone to take a while
  prio <- prioritze.by.randomwalk(locus_graph, 
                          joined_graph_ppi, 
                          nodes(ppi_graph),
                          tf_nodes)
  
  if(gene_only) {
    prio$best.trans.gene
  } else {
    prio
  }
}

