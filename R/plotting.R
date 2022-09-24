#' Plot a network
#'
#' Plot a network using degree-modulated node sizes, community colors and other enhancements
#' @param gr the network as an igraph object
#' @param community community assignment; vector of node labels
#' @param ... other settings
#' @export
plot_net = function(gr,
                    community = NULL,
                    color_map = NULL,
                    #deg_perc = 0.997,
                    extract_lcc = TRUE,
                    heavy_edge_deg_perc = 0.97,
                    coord = NULL,
                    vsize_func = function(deg) log(deg+3)*1, vertex_border = F,
                    niter = 1000, # number of iteration for FR layout computation
                    vertex_alpha = 0.4,
                    remove_loops = T,
                    make_simple = F,
                    show_legend = F, ...) {

  check_pkg_and_stop("igraph", "plot_net")

  gr = igraph::simplify(gr, remove.multiple = F, remove.loops = TRUE)
  igraph::V(gr)$label = NA

  if (is.null(community)) {
    igraph::V(gr)$color = adjustcolor("red", vertex_alpha)
  } else {
    if (is.null(color_map)) {
      # original color
      color_map = rainbow(max(community))
      # new color
      # color_map <- c("#56B4E9", "#CC79A7", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00")

    }
    igraph::V(gr)$color =  adjustcolor(color_map[community], vertex_alpha)
  }

  if (extract_lcc) {
    # Extract largest connected component after setting community colors to avoid losing community information
    gr = extract_largest_cc(gr)
  }
  degs = igraph::degree(gr)

  heavy_edge_thresh = quantile(degs, heavy_edge_deg_perc)
  e_wei = apply(igraph::as_edgelist(gr, names=F), 1, function(x) min(degs[x[1]], degs[x[2]]))
  eidx = e_wei > heavy_edge_thresh

  #arrow_size = rep(0.1, ecount(gr))
  if (make_simple) {
    igraph::E(gr)$weight = 1
    gr = igraph::simplify(gr, edge.attr.comb=list(weight="sum"))
    max_edge_wei = max(igraph::E(gr)$weight)
    igraph::E(gr)$width = 0.2*(1.3^igraph::E(gr)$weight)
    igraph::E(gr)$color = sapply(igraph::E(gr)$weight/max_edge_wei, function(x) adjustcolor("black", x))

  } else {
    igraph::E(gr)$width = 0.5
    igraph::E(gr)[eidx]$width = 1
    igraph::E(gr)$color = adjustcolor("gray",0.4)
    if (sum(eidx) > 0) {
      te_wei =  sqrt(e_wei[eidx]) # log(e_wei[eidx])
      igraph::E(gr)[eidx]$color = sapply(te_wei/max(te_wei), function(x) adjustcolor("darkblue", x)) # "darkblue"
    }
  }



  if (is.null(coord)) {
    coord <- igraph::layout_with_fr(gr, niter = niter)
  }

  plot(gr,
       layout = coord,
       #vertex.color = vcolors,
       #vertex.label = NA,
       vertex.size = vsize_func(degs),
       vertex.label.font = 1,
       vertex.label.cex = .7,
       vertex.label.color = "black",
       vertex.frame.color = ifelse(vertex_border, "black", NA),
       edge.arrow.size = 0.2) #, ... )

  if (show_legend)
    legend('topleft', legend=paste(1:N, ":", igraph::V(gr)[idx]$name),
           pt.cex=0.5*(log(degs[idx])), cex = 0.6, col=colors, pch=20, bty="n")

  list(gr = gr, coord = coord)
}


#' Plot degree distribution
#'
#' Plot the degree distribution of a network on log scale
#'
#' @param gr the network as an igraph object
#' @export
plot_deg_dist = function(gr, logx = T) {
  check_pkg_and_stop("igraph", "plot_deg_dist")

  deg_dist = igraph::degree_distribution(gr) # first element is the density of 0-degree nodes
  max_deg = length(deg_dist) - 1

  if (logx) {
    if (deg_dist[1] > 0) warning("There are 0-degree nodes. Omitting them on log scale.")
    x = 1:max_deg
    y = deg_dist[-1]
    xlabel = "Degree (log scale)"
    trans = "log10"
  } else {
    x = 0:max_deg
    y = deg_dist
    xlabel = "Degree"
    trans = "identity"
  }

  ggplot2::ggplot(data.frame(x = x, y = y)) +
    ggplot2::geom_segment(ggplot2::aes(x, y, xend=x, yend=0), color="slateblue") +
    ggplot2::scale_x_continuous(trans=trans) +
    ggplot2::labs(x=xlabel, y="Density", title="Degree Distribution") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank())
    # ggplot2::scale_y_continuous(expand=c(0,0), trans="sqrt") +
    # ggplot2::scale_x_continuous(trans='log10') +
    # ggplot2::labs(x="Degree (log scale)", y="Density (sqrt scale)", title="Degree Distribution") +
}

# hist(theta)
# lambda = 20
# b =  lambda / (n*(mean(theta)^2))
# dd = rpois(n, b*theta*sum(theta))
# summary(dd)
# max(dd)
# (hcount = hist(dd, seq(-0.5, max(dd)+0.5, by=1))$counts)
# ggplot2::ggplot(data.frame(x = 0:max(dd), y = hcount)) +
#   ggplot2::geom_segment(ggplot2::aes(x, y, xend=x, yend=0), color="slateblue") +
#   ggplot2::scale_x_continuous(trans="identity") +
#   ggplot2::labs(x="Degree", y="Density", title="Degree Distribution") +
#   ggplot2::theme_bw() +
#   ggplot2::theme(panel.border = ggplot2::element_blank())
#
# summary(dd)


