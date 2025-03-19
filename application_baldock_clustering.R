library("colSBM")
library("future")
library("future.apply")
library("future.callr")
library("aricode")
library("here")
options(future.globals.maxSize = Inf)

save_path <- here("applications", "baldock")
if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
}
start_time <- format(Sys.time(), "%Y%m%d%H%M%S")
temp_path <- file.path(save_path, paste0("tmp", start_time))
if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
}

all_dore_matrices <- readRDS(here("data", "dore-binary-matrices.Rds"))

# Filter to keep only Baldock

baldock_dore_matrices <- all_dore_matrices[grepl("Baldock", x = names(all_dore_matrices))]

conditions <- expand.grid(
    model = c("iid", "pi", "rho", "pirho"),
    rep = seq_len(5)
)
conditions$model <- as.character(conditions$model)
set.seed(1234)
plan(list(
    tweak("callr", workers = floor(parallelly::availableCores(omit = 1L) / 3L)),
    tweak("callr", workers = 3L)
))
clustering_results <- future_lapply(seq_len(nrow(conditions)), function(s) {
    model <- conditions[s, "model"]
    rep <- conditions[s, "rep"]
    message(
        "Starting simulation ", s, " on ", nrow(conditions),
        " with model = ", model, " and rep = ", rep
    )

    out <- clusterize_bipartite_networks(
        netlist = baldock_dore_matrices,
        colsbm_model = model,
        global_opts = list(backend = "future"),
        fit_opts = list(max_vem_steps = 5000L)
    )

    message("End simulation ", s, " on ", nrow(conditions))
    saveRDS(out, file = file.path(temp_path, paste0("simulation_", s, "_on_", nrow(conditions), "_", model, "_", rep, ".Rds")))
}, future.seed = TRUE)

saveRDS(clustering_results, file = file.path(save_path, paste0("clustering_results_", start_time, ".Rds")))
