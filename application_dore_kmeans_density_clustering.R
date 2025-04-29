library("colSBM")
library("future")
library("future.apply")
library("future.callr")
library("aricode")
library("here")
options(future.globals.maxSize = Inf)

save_path <- here("applications", "dore-kmeans-density")
if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
}
start_time <- format(Sys.time(), "%Y%m%d%H%M%S")
temp_path <- file.path(save_path, paste0("tmp", start_time))
if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
}

all_dore_matrices <- readRDS(here("data", "dore-binary-matrices.Rds"))

dore_densities <- sapply(all_dore_matrices, function(mat) {
    mat <- as.matrix(mat)
    density <- sum(mat) / (nrow(mat) * ncol(mat))
    return(density)
})

library(factoextra)

# Plot the silhouette method to find the optimal number of clusters
fviz_nbclust(as.matrix(dore_densities), kmeans, method = "silhouette")

precluster_density <- kmeans(dore_densities, centers = 2, nstart = 25)$cluster

first_matrices <- all_dore_matrices[precluster_density == 1]
second_matrices <- all_dore_matrices[precluster_density == 2]


conditions <- expand.grid(
    model = c("iid", "pi", "rho", "pirho"),
    rep = seq_len(5)
)
conditions$model <- as.character(conditions$model)
set.seed(1234)
plan(list(
    tweak("callr", workers = floor((parallelly::availableCores(omit = 1L) / 6L) / 2L)),
    tweak("callr", workers = 2L),
    tweak("callr", workers = 3L)
))
clustering_results <- future_lapply(seq_len(nrow(conditions)), function(s) {
    model <- conditions[s, "model"]
    rep <- conditions[s, "rep"]
    message(
        "Starting condition ", s, " on ", nrow(conditions),
        " with model = ", model, " and rep = ", rep
    )

    prefit <- c(
        estimate_colBiSBM(
            netlist = first_matrices,
            colsbm_model = model,
            global_opts = list(backend = "future"),
            fit_opts = list(max_vem_steps = 5000L)
        ),
        estimate_colBiSBM(
            netlist = second_matrices,
            colsbm_model = model,
            global_opts = list(backend = "future"),
            fit_opts = list(max_vem_steps = 5000L)
        )
    )

    out <- clusterize_bipartite_networks(
        netlist = all_dore_matrices,
        colsbm_model = model,
        global_opts = list(backend = "future"),
        fit_opts = list(max_vem_steps = 5000L),
        partition_init = prefit
    )

    message("End simulation ", s, " on ", nrow(conditions))
    saveRDS(out, file = file.path(temp_path, paste0("simulation_", s, "_on_", nrow(conditions), "_", model, "_", rep, ".Rds")))
}, future.seed = TRUE)

saveRDS(clustering_results, file = file.path(save_path, paste0("clustering_results_", start_time, ".Rds")))
