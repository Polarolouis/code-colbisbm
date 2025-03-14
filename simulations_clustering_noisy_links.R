library(colSBM)
library(future)
library(future.apply)
library(future.callr)
library(aricode)
library(here)
options(future.globals.maxSize = Inf)

save_path <- here("simulations", "clustering", "noisy_links")
if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
}

start_time <- format(Sys.time(), "%Y%m%d%H%M%S")

temp_path <- here("simulations", "clustering", "noisy_links", paste0("tmp", start_time))
if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
}
nr <- 120
nc <- 120

pi1 <- c(0.5, 0.3, 0.2)
rho1 <- c(0.4, 0.3, 0.2, 0.1)
alpha1 <- matrix(c(
    0.85, 0.4, 0.2, 0.05,
    0.6, 0.2, 0.05, 0.05,
    0.2, 0.05, 0.05, 0.7
), nrow = 3L)

pi2 <- c(0.5, 0.3, 0.2)
rho2 <- c(0.45, 0.3, 0.25)
alpha2 <- matrix(c(
    0.65, 0.05, 0.05,
    0.05, 0.8, 0.05,
    0.05, 0.05, 0.4
), nrow = 3L)

M <- 15L

eps_seq <- seq(0, 0.5, 0.05)
n_reps <- 5L
rep_clustering <- 5L

conditions <- expand.grid(
    epsilon = eps_seq,
    rep = seq_len(n_reps)
)

plan(list(tweak("callr", workers = floor(parallelly::availableCores(omit = 1L) / rep_clustering * 3L)), tweak("callr", workers = 5L), tweak("callr", workers = 3L)))

set.seed(1234)
df_list <- future_lapply(seq_len(nrow(conditions)), function(s) {
    eps <- conditions[s, "epsilon"]
    rep <- conditions[s, "rep"]

    coll_type1 <- generate_bipartite_collection(
        nr = nr,
        nc = nc,
        pi = pi1,
        rho = rho1,
        alpha = alpha1,
        M = M
    )
    names(coll_type1) <- paste(rep("type1", M), seq(1, M), sep = ".")

    coll_type2 <- generate_bipartite_collection(
        nr = nr,
        nc = nc,
        pi = pi2,
        rho = rho2,
        alpha = alpha2,
        M = M
    )
    names(coll_type2) <- paste(rep("type2", M), seq(1, M), sep = ".")

    coll <- c(coll_type1, coll_type2)

    noisy_coll <- lapply(coll, function(mat) {
        links <- sample.int(nr * nc, size = floor(eps * nr * nc))
        mat[links] <- as.integer(!mat[links])
        mat
    })

    type_clust <- rep(c("noisy", "clear"), rep_clustering)

    clustering_results <- future_lapply(type_clust, function(type) {
        current_netlist <- switch(type,
            "noisy" = noisy_coll,
            "clear" = coll,
            stop("Invalid type")
        )
        clust <- clusterize_bipartite_networks(
            netlist = current_netlist,
            colsbm_model = "iid",
            net_id = names(current_netlist),
            global_opts = list(backend = "future"),
            fit_opts = list(max_vem_steps = 5000L)
        )
        clust[c("cluster", "elapsed_time")]
    }, future.seed = TRUE)

    # TODO Add partial saving of results

    # Compute ARI to ground truth
    ground_truth <- rep(c(1, 2), each = M)
    ari_results <- sapply(clustering_results, function(res) {
        ari <- ARI(ground_truth, res$cluster)
    })
    elapsed_time <- sapply(clustering_results, function(res) {
        units(res$elapsed_time) <- "secs"
        res$elapsed_time
    })

    current_dataset_df <- data.frame(
        epsilon = eps,
        condition_rep = rep,
        clustering_rep = seq(1, rep_clustering),
        type = type_clust,
        ari = ari_results,
        elapsed_time = elapsed_time
    )
    saveRDS(current_dataset_df, file = file.path(temp_path, paste0("simulations_clustering_noisy_links_", start_time, "_", s, ".rds")))
}, future.seed = TRUE)

df <- do.call(rbind, df_list)
saveRDS(df, file = file.path(save_path, paste0("simulations_clustering_noisy_links_", start_time, ".rds")))
