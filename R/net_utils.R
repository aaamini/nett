# Extract a low degree connected component of a graph
#' @export
extract_low_deg_comp = function(g, deg_prec = 0.75, verb =F) {
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
           vcount(g2),
           vcount(g)/vcount(g2),
           round(max(degs) / max(degs2)),
           round(mean(degs) / mean(degs2)))
  }

  return(g2)
}


# Extract the largest connected component of a graph
#' @export
extract_largest_cc = function(gr, mode = "weak") {
  cc = igraph::components(gr, mode = mode)
  max_cc_idx = which.max(cc$csize)
  igraph::induced_subgraph(gr, cc$membership == max_cc_idx)
}
