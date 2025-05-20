library(colSBM)
library(here)
library(stringr)
library(future.apply)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 8912896000) # 8.3 GB

baldock_matrices <- readRDS(here("data", "baldock2019-binary-matrices.Rds"))

set.seed(123L)
if (!file.exists(here("data", "baldock-iid.Rds"))) {
    fit <- estimate_colBiSBM(
        netlist = baldock_matrices,
        net_id = names(baldock_matrices),
        colsbm_model = "iid",
        global_opts = list(
            verbosity = 0,
            backend = "future"
        ),
        fit_opts = list(
            max_vem_steps = 10000L
        )
    )
    saveRDS(fit, here("data", "baldock-iid.Rds"))
} else {
    fit <- readRDS(here("data", "baldock-iid.Rds"))
}

# for (i in seq_along(fit$sep_BiSBM$models)) {
#     cat("Network ID:", fit$sep_BiSBM$models[[i]]$net_id, "\n")
#     pdf(
#         here("figures", "applications", "baldock", paste0(
#             "bisbm-mat-",
#             fit$sep_BiSBM$models[[i]]$net_id,
#             ".pdf"
#         )),
#         family = "Times"
#     )
#     print(plot(fit$sep_BiSBM$models[[i]], type = "block") +
#         xlab("Plants") +
#         ylab("Pollinators") +
#         theme(axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40)))
#     dev.off()
# }

# for (i in seq_along(fit$sep_BiSBM$models)) {
#     cat("Network ID:", fit$sep_BiSBM$models[[i]]$net_id, "\n")
#     pdf(
#         here("figures", "applications", "baldock", paste0(
#             "mat-",
#             fit$sep_BiSBM$models[[i]]$net_id,
#             ".pdf"
#         )),
#         family = "Times"
#     )
#     print(sbm::plotMyMatrix(fit$sep_BiSBM$models[[i]]$A[[1]]) +
#         xlab("Plants") +
#         ylab("Pollinators") +
#         guides(fill = "none") +
#         theme_bw() +
#         coord_fixed(ratio = 1) +
#         theme(axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40), axis.text = element_blank(), axis.ticks = element_blank(), strip.text = element_blank(), strip.background = element_rect(fill = "white", linewidth = 0)))
#     dev.off()
# }

# for (i in seq_along(fit$sep_BiSBM$models)) {
#     net_id <- fit$sep_BiSBM$models[[i]]$net_id
#     cat("Network ID:", net_id, "\n")
#     print(plot(fit$sep_BiSBM$models[[i]], type = "block") +
#         xlab("Plants") +
#         ylab("Pollinators") +
#         theme(axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40)))
# }
# plot(fit)

# plot(fit$best_fit, type = "meso", mixture = TRUE, values = TRUE)

baldock_nodes_groups <- extract_nodes_groups(fit)
baldock_nodes_groups <- baldock_nodes_groups %>% mutate(node_type = ifelse(node_type == "row", "pollinator", "plant"), network = recode(network,
    "1" = "Bristol",
    "2" = "Edinburgh",
    "3" = "Leeds",
    "4" = "Reading"
))

netnames <- c("Bristol", "Edinburgh", "Leeds", "Reading")

# Investigate species present in the 4 networks

shared_nodes <- baldock_nodes_groups %>%
    group_by(node_name, node_type) %>%
    summarise(n_networks = n_distinct(network), .groups = "drop") %>%
    filter(n_networks == 4)

baldock_shared_nodes <- baldock_nodes_groups %>%
    semi_join(shared_nodes, by = c("node_name", "node_type"))

nodes_diff_cluster <- baldock_shared_nodes %>%
    group_by(node_name, node_type) %>%
    summarise(n_clusters = n_distinct(cluster), .groups = "drop") %>%
    filter(n_clusters > 1)
baldock_shared_nodes_diff_cluster <- baldock_shared_nodes %>%
    semi_join(nodes_diff_cluster, by = c("node_name", "node_type"))

baldock_shared_nodes_diff_cluster %>%
    group_by(node_name, node_type) %>%
    arrange(.by_group = TRUE) %>%
    View()

bombus_1 <- which(rownames(baldock_matrices[[1]]) |> grepl(pattern = "lapidarius|hortorum"))
bombus_2 <- which(rownames(baldock_matrices[[2]]) |> grepl(pattern = "lapidarius|hortorum"))
bombus_3 <- which(rownames(baldock_matrices[[3]]) |> grepl(pattern = "lapidarius|hortorum"))
bombus_4 <- which(rownames(baldock_matrices[[4]]) |> grepl(pattern = "lapidarius|hortorum"))

bombus_matrix_bristol <- baldock_matrices[[1]][bombus_1, ] |> bipartite::empty()
baldock_matrices[[2]][bombus_2, ] |> bipartite::empty() -> bombus_matrix_edinburgh
baldock_matrices[[3]][bombus_3, ] |> bipartite::empty() -> bombus_matrix_leeds
baldock_matrices[[4]][bombus_4, ] |> bipartite::empty() -> bombus_matrix_reading

bombus_mat_list <- list(
    bombus_matrix_bristol,
    bombus_matrix_edinburgh,
    bombus_matrix_leeds,
    bombus_matrix_reading
)

library(tidygraph)
library(ggraph)
library(ggokabeito)
library(patchwork)
lapply(seq_along(bombus_mat_list), function(idx) {
    bombus_mat_list[[idx]] |>
        as.matrix() |>
        as_tbl_graph() |>
        mutate(node_type = ifelse(type, "Plant", "Pollinator"), bombus = ifelse(node_type == "Pollinator", ifelse(grepl(pattern = "hortorum", name), "Hortorum", "Lapidarius"), "Plant")) |>
        ggraph(layout = "igraph", algorithm = "bipartite") +
        geom_node_point(aes(color = node_type, shape = bombus), size = 10) +
        labs(shape = "Species") +
        geom_node_text(aes(label = ifelse(bombus != "Plant", bombus, "")), nudge_y = 0.05) +
        geom_edge_link(alpha = 0.1) +
        scale_color_okabe_ito(order = 6:7) +
        theme_void() +
        ggtitle(netnames[idx]) +
        labs(color = "Node type") +
        theme(legend.position = "bottom")
}) |> wrap_plots(ncol = 2) +
    plot_annotation(title = "Bombus species interactions with plants in the four networks") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
