library(here)
library(ggplot2)
library(colSBM)


bipartite_graph_density <- function(A) {
    # Calculate the density of a bipartite graph
    # A: adjacency matrix of the bipartite graph
    # Returns the density of the bipartite graph

    n <- nrow(A)
    m <- ncol(A)
    d <- sum(A) / (n * m)
    return(d)
}

extract_vbound_per_net <- function(fit) {
    sapply(seq_along(fit$A), function(m) {
        fit$entropy_tau(m) + fit$vb_tau_pi(m) + fit$vb_tau_alpha(m)
    })
}

all_dore_matrices <- readRDS(here("data", "dore-binary-matrices.Rds"))
density_all_dore <- sapply(
    all_dore_matrices,
    bipartite_graph_density
)
density_all_dore_df <- data.frame(
    density = density_all_dore,
    name = names(density_all_dore)
)

# Filter to keep only Baldock Traveset Souza Cordeniz Trojelsgaard and Gibson

subdore_matrices <- all_dore_matrices[grepl("Baldock|Traveset|Souza|Cordeniz|Trojelsgaard|Gibson", x = names(all_dore_matrices))]

baldock_matrices <- all_dore_matrices[grepl("Baldock", x = names(all_dore_matrices))]

density_subdore <- sapply(
    subdore_matrices,
    bipartite_graph_density
)
library(tidyverse)
density_subdore_df <- data.frame(
    density = density_subdore,
    name = names(density_subdore)
)

ggplot(density_subdore_df) +
    aes(x = reorder(name, -density), y = density) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_x_discrete(labels = abbreviate) +
    labs(
        x = "Network",
        y = "Density",
        title = "Density of bipartite graphs for selected Doré networks"
    ) +
    theme_minimal()

ggplot(density_all_dore_df) +
    aes(x = reorder(name, -density), y = density) +
    geom_bar(stat = "identity") +
    scale_x_discrete(labels = abbreviate) +
    scale_y_continuous(n.breaks = 10) +
    coord_flip() +
    labs(
        x = "Network",
        y = "Density",
        title = "Density of bipartite graphs for all Doré networks"
    ) +
    theme_minimal()

density_subdore_df

# Trying to fit on baldock and souza density < 0.05

set.seed(123)
networks_density_0.05 <- subdore_matrices[density_subdore_df$density < 0.05]
library(sbm)

sep_fit <- lapply(networks_density_0.05, function(mat) {
    estimateBipartiteSBM(netMat = mat)
})

sep_BICL <- sum(sapply(sep_fit, function(x) {
    x$loglik + x$entropy - x$penalty
}))

sep_vbound_0.05 <- sapply(sep_fit, function(x) {
    x$loglik + x$entropy
})

fit_baldock <- clusterize_bipartite_networks(
    netlist = baldock_matrices,
    colsbm_model = "iid",
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "no_mc"
    )
)

fit_0.05 <- clusterize_bipartite_networks(
    netlist = networks_density_0.05,
    colsbm_model = "iid",
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "no_mc"
    )
)
vbound_0.05 <- extract_vbound_per_net(fit_0.05$partition[[1]])
ggplot(data.frame(network = names(sep_vbound_0.05), vb_ratio = vbound_0.05 / sep_vbound_0.05)) +
    aes(y = vb_ratio, x = network) +
    geom_bar(stat = "identity") +
    # coord_flip() +
    # scale_x_discrete(labels = abbreviate) +
    labs(
        x = "Network",
        y = "VB ratio",
        title = "VB ratio for networks with density < 0.05"
    ) +
    theme_minimal()


# density between 0.05 and 0.1
networks_density_0.1 <- subdore_matrices[density_subdore_df$density >= 0.05 & density_subdore_df$density < 0.1]

sep_fit_0.1 <- lapply(networks_density_0.1, function(mat) {
    estimateBipartiteSBM(netMat = mat)
})
sep_BICL_0.1 <- sum(sapply(sep_fit_0.1, function(x) {
    x$loglik + x$entropy - x$penalty
}))
sep_vbound_0.1 <- sapply(sep_fit_0.1, function(x) {
    x$loglik + x$entropy
})
fit_0.1 <- clusterize_bipartite_networks(
    netlist = networks_density_0.1,
    colsbm_model = "iid",
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "no_mc"
    )
)
vbound_0.1 <- extract_vbound_per_net(fit_0.1$partition[[1]])
ggplot(data.frame(network = names(sep_vbound_0.1), vb_ratio = vbound_0.1 / sep_vbound_0.1)) +
    aes(y = vb_ratio, x = network) +
    geom_bar(stat = "identity") +
    # coord_flip() +
    # scale_x_discrete(labels = abbreviate) +
    labs(
        x = "Network",
        y = "VB ratio",
        title = "VB ratio for networks with density < 0.1"
    ) +
    theme_minimal()

# density between 0.1 and 0.2
networks_density_0.2 <- subdore_matrices[density_subdore_df$density >= 0.1 & density_subdore_df$density < 0.2]
sep_fit_0.2 <- lapply(networks_density_0.2, function(mat) {
    estimateBipartiteSBM(netMat = mat)
})
sep_BICL_0.2 <- sum(sapply(sep_fit_0.2, function(x) {
    x$loglik + x$entropy - x$penalty
}))
sep_vbound_0.2 <- sapply(sep_fit_0.2, function(x) {
    x$loglik + x$entropy
})
fit_0.2 <- clusterize_bipartite_networks(
    netlist = networks_density_0.2,
    colsbm_model = "iid",
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "no_mc"
    )
)
vbound_0.2 <- extract_vbound_per_net(fit_0.2$partition[[1]])
ggplot(data.frame(network = names(sep_vbound_0.2), vb_ratio = vbound_0.2 / sep_vbound_0.2)) +
    aes(y = vb_ratio, x = network) +
    geom_bar(stat = "identity") +
    # coord_flip() +
    # scale_x_discrete(labels = abbreviate) +
    labs(
        x = "Network",
        y = "VB ratio",
        title = "VB ratio for networks with density < 0.2"
    ) +
    theme_minimal()


# density between 0.2 and 0.3
networks_density_0.3 <- subdore_matrices[density_subdore_df$density >= 0.2 & density_subdore_df$density < 0.3]
sep_fit_0.3 <- lapply(networks_density_0.3, function(mat) {
    estimateBipartiteSBM(netMat = mat)
})
sep_BICL_0.3 <- sum(sapply(sep_fit_0.3, function(x) {
    x$loglik + x$entropy - x$penalty
}))
sep_vbound_0.3 <- sapply(sep_fit_0.3, function(x) {
    x$loglik + x$entropy
})
fit_0.3 <- clusterize_bipartite_networks(
    netlist = networks_density_0.3,
    colsbm_model = "iid",
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "no_mc"
    )
)
vbound_0.3 <- c(extract_vbound_per_net(fit_0.3$partition[[1]]), extract_vbound_per_net(fit_0.3$partition[[2]]))
ggplot(data.frame(network = names(sep_vbound_0.3), vb_ratio = vbound_0.3 / sep_vbound_0.3)) +
    aes(y = vb_ratio, x = network) +
    geom_bar(stat = "identity") +
    # coord_flip() +
    # scale_x_discrete(labels = abbreviate) +
    labs(
        x = "Network",
        y = "VB ratio",
        title = "VB ratio for networks with density < 0.3"
    ) +
    theme_minimal()

# Density between 0.05 and 0.1
networks_density_0.05p0.01 <- c(networks_density_0.05, networks_density_0.1)
sep_fit_0.05p0.01 <- lapply(networks_density_0.05p0.01, function(mat) {
    estimateBipartiteSBM(netMat = mat)
})
sep_BICL_0.05p0.01 <- sum(sapply(sep_fit_0.05p0.01, function(x) {
    x$loglik + x$entropy - x$penalty
}))
sep_vbound_0.05p0.01 <- sapply(sep_fit_0.05p0.01, function(x) {
    x$loglik + x$entropy
})

vbound_0.05p0.01 <- c(extract_vbound_per_net(fit_0.05p0.01$partition[[1]]), extract_vbound_per_net(fit_0.05p0.01$partition[[2]]))
ggplot(data.frame(network = names(sep_vbound_0.05p0.01), vb_ratio = vbound_0.05p0.01 / sep_vbound_0.05p0.01)) +
    aes(y = vb_ratio, x = network) +
    geom_bar(stat = "identity") +
    # coord_flip() +
    # scale_x_discrete(labels = abbreviate) +
    labs(
        x = "Network",
        y = "VB ratio",
        title = "VB ratio for networks with density < 0.1"
    ) +
    theme_minimal()

library(future.apply)
library(future.callr)
options(future.globals.maxSize = Inf)
plan(tweak("callr", workers = 10L))
clust_res_0.05p0.01 <- future_sapply(seq(1, 10), function(idx) {
    fit_0.05p0.01 <- clusterize_bipartite_networks(
        netlist = networks_density_0.05p0.01,
        colsbm_model = "iid",
        distribution = "bernoulli",
        fit_opts = list(
            algo_ve = "fp",
            minibatch = TRUE,
            verbosity = 0
        ),
        global_opts = list(
            verbosity = 2,
            backend = "no_mc"
        )
    )
    fit_0.05p0.01$cluster
}, future.seed = TRUE)


#
clust_res_0.05 <- future_sapply(seq(1, 10), function(idx) {
    fit <- clusterize_bipartite_networks(
        netlist = networks_density_0.05,
        colsbm_model = "iid",
        distribution = "bernoulli",
        fit_opts = list(
            algo_ve = "fp",
            minibatch = TRUE,
            verbosity = 0
        ),
        global_opts = list(
            verbosity = 2,
            backend = "no_mc"
        )
    )
    fit$cluster
})
clust_res_0.1 <- future_sapply(seq(1, 10), function(idx) {
    fit <- clusterize_bipartite_networks(
        netlist = networks_density_0.1,
        colsbm_model = "iid",
        distribution = "bernoulli",
        fit_opts = list(
            algo_ve = "fp",
            minibatch = TRUE,
            verbosity = 0
        ),
        global_opts = list(
            verbosity = 2,
            backend = "no_mc"
        )
    )
    fit$cluster
}, future.seed = TRUE)
clust_res_0.2 <- future_sapply(seq(1, 10), function(idx) {
    fit <- clusterize_bipartite_networks(
        netlist = networks_density_0.2,
        colsbm_model = "iid",
        distribution = "bernoulli",
        fit_opts = list(
            algo_ve = "fp",
            minibatch = TRUE,
            verbosity = 0
        ),
        global_opts = list(
            verbosity = 2,
            backend = "no_mc"
        )
    )
    fit$cluster
}, future.seed = TRUE)
clust_res_0.3 <- future_sapply(seq(1, 10), function(idx) {
    fit <- clusterize_bipartite_networks(
        netlist = networks_density_0.3,
        colsbm_model = "iid",
        distribution = "bernoulli",
        fit_opts = list(
            algo_ve = "fp",
            minibatch = TRUE,
            verbosity = 0
        ),
        global_opts = list(
            verbosity = 2,
            backend = "no_mc"
        )
    )
    fit$cluster
}, future.seed = TRUE)


fit_subdore <- colSBM::estimate_colBiSBM(
    netlist = subdore_matrices,
    net_id = names(subdore_matrices),
    colsbm_model = "iid",
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "no_mc"
    )
)

vbound_subdore <- extract_vbound_per_net(fit_subdore$best_fit)
sep_vbound_subdore <- sapply(fit_subdore$sep_BiSBM$models, function(x) x$vb_tau_pi(1) + x$vb_tau_alpha(1) + x$entropy_tau(1))
ggplot(data.frame(network = names(subdore_matrices), vb_ratio = vbound_subdore / sep_vbound_subdore)) +
    aes(y = vb_ratio, x = network) +
    geom_bar(stat = "identity") +
    # coord_flip() +
    # scale_x_discrete(labels = abbreviate) +
    labs(
        x = "Network",
        y = "VB ratio",
        title = "VB ratio for networks subdore"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
fit <- fit_subdore
dist_bm <- compute_dissimilarity_matrix(collection = fit)
cl_kmed <- partition_networks_list_from_dissimilarity(
    networks_list = fit$A,
    dissimilarity_matrix = dist_bm,
    nb_groups = 2L
)
dist_bm <- sqrt(dist_bm)
hc_wardD2 <- hclust(as.dist(dist_bm), method = "ward.D2")
cl_wardD2 <- cutree(hc_wardD2, k = 2L)
hc_ward <- hclust(as.dist(dist_bm), method = "ward.D")
cl_ward <- cutree(hc_ward, k = 2L)
hc_single <- hclust(as.dist(dist_bm), method = "single")
(cl_single <- cutree(hc_single, k = 2L))
hc_complete <- hclust(as.dist(dist_bm), method = "complete")
(cl_complete <- cutree(hc_complete, k = 2L))
hc_average <- hclust(as.dist(dist_bm), method = "average")
(cl_average <- cutree(hc_average, k = 2L))
hc_centroid <- hclust(as.dist(dist_bm), method = "centroid")
(cl_centroid <- cutree(hc_centroid, k = 2L))
hc_mcquitty <- hclust(as.dist(dist_bm), method = "mcquitty")
plot(hc_mcquitty)
(cl_mcquitty <- cutree(hc_mcquitty, k = 2L))
hc_median <- hclust(as.dist(dist_bm), method = "median")
plot(hc_median)
(cl_median <- cutree(hc_median, k = 2L))

hcdiana <- cluster::diana(as.dist(dist_bm))

unname(cl_kmed)

cluster_subdore_rep <- future_lapply(seq(1, 10), function(idx) {
    fit <- clusterize_bipartite_networks(
        netlist = subdore_matrices,
        colsbm_model = "iid",
        net_id = names(subdore_matrices),
        distribution = "bernoulli",
        fit_opts = list(
            algo_ve = "fp",
            minibatch = TRUE,
            verbosity = 0
        ),
        global_opts = list(
            verbosity = 0,
            backend = "no_mc"
        )
    )
    fit
}, future.seed = TRUE)

saveRDS(cluster_subdore_rep, file = here("data", "cluster_subdore_rep.Rds"))
