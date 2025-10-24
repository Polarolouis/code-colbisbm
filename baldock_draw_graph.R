library("ggraph")
library("tidyr")
library("tidygraph")
library("tikzDevice")


middle_red_hex <- "#d55d52"
middle_blue_hex <- "#4894c4"

baldock_matrices_list <- readRDS("data/baldock2019-binary-matrices.Rds")
baldock_graphs_list <- lapply(baldock_matrices_list, function(mat) {
    as_tbl_graph(mat) |> mutate(node_type = ifelse(type, "Insect", "Plant"), type = !type)
})

graph_plot_list <- lapply(baldock_graphs_list, function(graph) {
    ggraph(graph = as_tbl_graph(graph), layout = "bipartite") +

        geom_edge_link(alpha = 0.4) +
        geom_node_point(aes(shape = node_type, color = node_type), size = 4) +
        scale_color_manual(values = c(middle_blue_hex, middle_red_hex)) +
        scale_shape_manual(values = c(16, 15)) +
        theme_void() +
        labs(color = "Node type", shape = "Node type") +
        theme(legend.position = "bottom", legend.direction = "horizontal")
})

base_tikz_path <- "tikz/applications/baldock"

options(tikzDocumentDeclaration = "\\documentclass{standalone}")
for (idx in seq_along(graph_plot_list)) {
    tikz(height = 3, file.path(base_tikz_path, paste0("graph-", names(graph_plot_list)[[idx]], ".tex")), standAlone = TRUE)
    print(graph_plot_list[[idx]])
    dev.off()
}
