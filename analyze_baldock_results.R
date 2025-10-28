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

for (i in seq_along(fit$sep_BiSBM$models)) {
    cat("Network ID:", fit$sep_BiSBM$models[[i]]$net_id, toString(fit$sep_BiSBM$models[[i]]$Q), "\n")
    pdf(
        here("figures", "applications", "baldock", paste0(
            "bisbm-mat-",
            fit$sep_BiSBM$models[[i]]$net_id,
            ".pdf"
        )),
        family = "Times"
    )
    print(plot(fit$sep_BiSBM$models[[i]], type = "block", oRow = seq(3, 1)) +
        xlab("Plants") +
        ylab("Pollinators") +
        theme(axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40)))
    dev.off()
}

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

for (i in seq_along(fit$net_id)) {
    net_id <- fit$net_id[i]
    cat("Network ID:", net_id, "\n")
    pdf(
        here("figures", "applications", "baldock", paste0(
            "colbisbm-mat-",
            fit$best_fit$net_id[i],
            ".pdf"
        )),
        family = "Times"
    )
    print(plot(fit$best_fit, net_id = i, type = "block", values = F) +
        xlab("Plants") +
        ylab("Pollinators") +
        theme(axis.title.y = element_text(size = 40), axis.title.x = element_text(size = 40)))
    dev.off()
}

# # plot(fit)
cairo_pdf(
    here("figures", "applications", "baldock", "shared-iid.pdf"),
    family = "Times", height = 5
)
plot(fit$best_fit, type = "meso", values = TRUE) + guides(fill = "none")
dev.off()
cairo_pdf(
    here("figures", "applications", "baldock", "shared-mixture-iid.pdf"),
    family = "Times", height = 7, width = 10
)
plot(fit$best_fit, type = "meso", mixture = TRUE, values = TRUE)
dev.off()
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

plants_bristol <- baldock_matrices[[1]] |> colnames()
plants_leeds <- baldock_matrices[[3]] |> colnames()
plants_b_and_l <- intersect(plants_bristol, plants_leeds)
length(plants_b_and_l)

plants_edinburgh <- baldock_matrices[[2]] |> colnames()
plants_reading <- baldock_matrices[[4]] |> colnames()
plants_e_and_r <- intersect(plants_edinburgh, plants_reading)
length(plants_e_and_r)

plants_all <- intersect(plants_b_and_l, plants_e_and_r)
length(plants_all)


design <-
    c(area(2, 2, 4, 4), area(1, 2, r = 4), area(2, 1, 4))
plot(design)


my_meso_plot <- function(best_fit) {
    plot(best_fit, type = "meso", values = TRUE, mixture = TRUE, text_size = 5.5) + plot_layout(guides = "collect") & theme(legend.position = "bottom", text = element_text(family = "serif"), legend.title = element_text(size = 12), axis.text = element_text(size = 11))
}

datalist <- readRDS(here("applications", "baldock", "baldock.Rds"))
fit_iid_1 <- datalist$iid.3$partition[[1]]
fit_iid_2 <- datalist$iid.3$partition[[2]]
class(fit_iid_1) <- c("fitBipartiteSBMPop", "R6")
class(fit_iid_2) <- c("fitBipartiteSBMPop", "R6")
fit_iid_1$net_id <- str_c("Kenya-", fit_iid_1$net_id |> str_replace_all("Baldock2011_", ""))
fit_iid_2$net_id <- fit_iid_2$net_id |> str_replace_all("Baldock201[1,9]_", "Brit-")
design_mat <- "AAABBBBB"
# (iid_plot <- wrap_plots("A" = my_meso_plot(fit_iid_1), "B" = my_meso_plot(fit_iid_2)) + plot_layout(guides = "collect", design = design_mat) &
#     theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical"))

(iid_plot <- my_meso_plot(fit_iid_2) & guides(fill = "none"))


ggsave("figures/applications/baldock/baldock-iid-clust.pdf", iid_plot)
fit_pirho_1 <- datalist$pirho.4$partition[[1]]

class(fit_pirho_1) <- c("fitBipartiteSBMPop", "R6")
fit_pirho_1$net_id <- fit_pirho_1$net_id |>
    str_replace("Baldock2011_", "Kenya-") |>
    str_replace("Baldock2011_", "") |>
    str_replace_all("Baldock2019_", "Brit-")
(pirho_plot <- my_meso_plot(fit_pirho_1) & guides(fill = "none"))
ggsave("figures/applications/baldock/baldock-pirho-clust.pdf", pirho_plot)

library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")
tikz(file = "figures/applications/baldock/baldock-iid-clust.tex", width = 10.3, height = 7.53, standAlone = TRUE)
iid_plot
dev.off()
tikz(file = "figures/applications/baldock/baldock-pirho-clust.tex", width = 10.3, height = 7.53, standAlone = TRUE)
pirho_plot
dev.off()
