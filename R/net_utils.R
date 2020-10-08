# Extract a low degree connected component of a graph
#' @export
extract_low_deg_comp = function(g, deg_prec = 0.75) {
  g = extract_largest_cc(g)
  degs = degree(g)
  g %>%
    induced_subgraph(degs <= quantile(degs, deg_prec))%>%
    extract_largest_cc() %>%
    as_adj
}


# Extract the largest connected component of a graph
#' @export
extract_largest_cc = function(gr, mode = "weak") {
  cc = igraph::components(gr, mode = mode)
  max_cc_idx = which.max(cc$csize)
  igraph::induced_subgraph(gr, cc$membership == max_cc_idx)
}
