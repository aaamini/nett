## code to prepare `polblogs` dataset goes here

polblogs = igraph::read_graph("data-raw/polblogs.gml", format = "gml")
polblogs = polblogs %>%
  igraph::set_vertex_attr("community", value = V(polblogs)$value+1)
polblogs = igraph::remove.vertex.attribute(polblogs, "value")
usethis::use_data(polblogs, overwrite = TRUE)
