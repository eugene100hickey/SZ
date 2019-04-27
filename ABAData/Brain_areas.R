setwd("~/Desktop/Academic/Research/Schiz/SZ")

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
require(ABAData)
library(ABAEnrichment)

data("dataset_adult")

area_numbers <- dataset_adult %>% distinct(structure)
area_names <- get_name(area_numbers$structure)
area_list <- map(area_numbers$structure, get_superstructures)
area_list_names <- map(area_list, get_name)

vertice1 <- as.character(unlist(area_list_names))
vertice2 <- vertice1
vertice1 <- vertice1[-length(vertice1)]
vertice2 <- vertice2[-1]
my_edges <- tibble(from = vertice1, to = vertice2)
my_edges <- my_edges %>% filter(to != "Br_Brain")
g <- graph_from_edgelist(as.matrix(my_edges), directed = T)


trim_and_plot <- function(graph = g, node = "Br_Brain", depth = 3){
  v_trimmed <- V(g)$name[distances(g, node, mode = "out") < depth]
  g1 <- as_tbl_graph(induced_subgraph(g, v_trimmed))
  text_colour <- distances(g1, node, mode = "out")
  text_colour <- as.numeric(text_colour)[order(as.numeric(text_colour))]
  ggraph(g1, 'dendrogram') + 
    geom_edge_link(colour = "lightblue") + 
    geom_node_text(aes(label = name, 
                       colour = text_colour, 
                        size = 8), angle = -10) + 
    xlim(-2, 10) + 
    ylim(-2, 3) +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
}
