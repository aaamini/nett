#' Extract low-degree component
#'
#' Extract a low-degree connected component of a network
#' @param g The network as an igraph object
#' @param deg_prec The cut-off degree percentile
#' @param verb Whether to be verbose (TRUE|FALSE)
#' @return An igraph object
#' @keywords utils
#' @export
extract_low_deg_comp = function(g, deg_prec = 0.75, verb =FALSE) {
  g = extract_largest_cc(g)
  degs = igraph::degree(g)
  g2 = g %>%
    igraph::induced_subgraph(degs <= quantile(degs, deg_prec)) %>%
    extract_largest_cc()

  if (verb) {
    degs2 = igraph::degree(g2)
    printf('new # of nodes = %d
reduction in # of nodes = %2.1fx
reduction in max deg. = %dx
reduction in mean deg. = %dx',
           igraph::vcount(g2),
           igraph::vcount(g)/igraph::vcount(g2),
           round(max(degs) / max(degs2)),
           round(mean(degs) / mean(degs2)))
  }

  return(g2)
}


#' Extract largest component
#'
#' Extract the largest connected component of a network
#' @param gr The network as an igraph object
#' @param mode Type of connected component ("weak"|"strong")
#' @return An igraph object
#' @keywords utils
#' @export
extract_largest_cc = function(gr, mode = "weak") {
  cc = igraph::components(gr, mode = mode)
  max_cc_idx = which.max(cc$csize)
  igraph::induced_subgraph(gr, cc$membership == max_cc_idx)
}
