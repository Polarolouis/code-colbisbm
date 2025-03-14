library(colSBM)
library(future)
library(future.apply)
library(future.callr)
library(aricode)
library(here)
options(future.globals.maxSize = Inf)

save_path <- here("simulations", "clustering", "noisy_alpha")
if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
}

start_time <- format(Sys.time(), "%Y%m%d%H%M%S")

temp_path <- here(save_path, paste0("tmp", start_time))
if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
}
nr <- 120
nc <- 120

pi1 <- c(0.5, 0.3, 0.2)
rho1 <- c(0.4, 0.3, 0.2, 0.1)
alpha1 <- matrix(c(
    0.85, 0.4, 0.2, 0.15,
    0.6, 0.2, 0.15, 0.15,
    0.2, 0.15, 0.15, 0.7
), nrow = 3L)

pi2 <- c(0.5, 0.3, 0.2)
rho2 <- c(0.45, 0.3, 0.25)
alpha2 <- matrix(c(
    0.65, 0.15, 0.15,
    0.15, 0.8, 0.15,
    0.15, 0.15, 0.4
), nrow = 3L)

M <- 15L

eps_seq <- seq(0, 0.05, 0.01)
n_reps <- 5L
rep_clustering <- 5L

conditions <- expand.grid(
    epsilon = eps_seq,
    rep = seq_len(n_reps)
)

plan(list(tweak("callr", workers = floor(parallelly::availableCores(omit = 1L) / (rep_clustering * 3L))), tweak("callr", workers = 5L), tweak("callr", workers = 3L)))

set.seed(1234)
df_list <- future_lapply(seq_len(nrow(conditions)), function(s) {
    eps <- conditions[s, "epsilon"]
    rep <- conditions[s, "rep"]
    message("Starting simulation ", s, " with epsilon = ", eps, " and rep = ", rep)
    # coll_type1 <- generate_bipartite_collection(
    #     nr = nr,
    #     nc = nc,
    #     pi = pi1,
    #     rho = rho1,
    #     alpha = alpha1,
    #     M = M
    # )

    coll_type1 <- lapply(seq(1, M), function(i) {
        noise_matrix <- matrix(rnorm(nrow(alpha1) * ncol(alpha1), sd = eps), nrow = nrow(alpha1))
        alpha1_noisy <- matrix(pmax(0, pmin(1, alpha1 + noise_matrix)), nrow = nrow(alpha1))
        generate_bipartite_collection(
            nr = nr,
            nc = nc,
            pi = pi1,
            rho = rho1,
            alpha = alpha1,
            M = 1
        )[[1]]
    })

    names(coll_type1) <- paste(rep("type1", M), seq(1, M), sep = ".")

    coll_type2 <- lapply(seq(1, M), function(i) {
        noise_matrix <- matrix(rnorm(nrow(alpha2) * ncol(alpha2), sd = eps), nrow = nrow(alpha2))
        alpha2_noisy <- matrix(pmax(0, pmin(1, alpha2 + noise_matrix)), nrow = nrow(alpha2))
        generate_bipartite_collection(
            nr = nr,
            nc = nc,
            pi = pi2,
            rho = rho2,
            alpha = alpha2,
            M = 1
        )[[1]]
    })
    names(coll_type2) <- paste(rep("type2", M), seq(1, M), sep = ".")

    coll <- c(coll_type1, coll_type2)

    clustering_results <- future_lapply(seq(1, rep_clustering), function(i) {
        clust <- clusterize_bipartite_networks(
            netlist = coll,
            colsbm_model = "iid",
            net_id = coll,
            global_opts = list(backend = "future"),
            fit_opts = list(max_vem_steps = 5000L)
        )
        clust[c("cluster", "elapsed_time")]
    }, future.seed = TRUE)

    # Compute ARI to ground truth
    ground_truth <- rep(c(1, 2), each = M)
    ari_results <- sapply(clustering_results, function(res) {
        ari <- ARI(ground_truth, res$cluster)
    })
    elapsed_time <- sapply(clustering_results, function(res) {
        units(res$elapsed_time) <- "secs"
        res$elapsed_time
    })

    ari_other_to_other <- outer(seq(1, rep_clustering), 1:rep_clustering, Vectorize(function(i, j) {
        ari <- ARI(clustering_results[[i]]$cluster, clustering_results[[j]]$cluster)
    }))
    colnames(ari_other_to_other) <- paste("ari_to",
        seq(1, rep_clustering),
        sep = "."
    )
    ari_other_to_other <- as.data.frame(ari_other_to_other)


    current_dataset_df <- data.frame(
        epsilon = eps,
        condition_rep = rep,
        clustering_rep = seq(1, rep_clustering),
        ari_truth = ari_results,
        elapsed_time = elapsed_time
    )
    current_dataset_df <- cbind(current_dataset_df, ari_other_to_other)
    saveRDS(current_dataset_df, file = file.path(temp_path, paste0("simulations_clustering_noisy_alpha_", start_time, "_", s, ".rds")))
    message("Finished simulation ", s, " with epsilon = ", eps, " and rep = ", rep)
}, future.seed = TRUE)

df <- do.call(rbind, df_list)
saveRDS(df, file = file.path(save_path, paste0("simulations_clustering_noisy_alpha_", start_time, ".rds")))
