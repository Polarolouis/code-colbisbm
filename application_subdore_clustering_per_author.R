library("colSBM")
library("future")
library("future.apply")
library("future.callr")
library("aricode")
library("here")
options(future.globals.maxSize = Inf)

save_path <- here("applications", "subdore")
if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
}
start_time <- format(Sys.time(), "%Y%m%d%H%M%S")
temp_path <- file.path(save_path, paste0("tmp", start_time))
if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
}

all_dore_matrices <- readRDS(here("data", "dore-binary-matrices.Rds"))

nb_rep <- 5L

# Filter to keep only Baldock Traveset Souza Cordeniz Trojelsgaard and Gibson

subdore_matrices <- all_dore_matrices[grepl("Baldock|Traveset|Souza|Cordeniz|Trojelsgaard|Gibson", x = names(all_dore_matrices))]

baldock_matrices <- all_dore_matrices[grepl("Baldock", x = names(all_dore_matrices))]
traveset_matrices <- all_dore_matrices[grepl("Traveset", x = names(all_dore_matrices))]
souza_matrices <- all_dore_matrices[grepl("Souza", x = names(all_dore_matrices))]
cordeniz_matrices <- all_dore_matrices[grepl("Cordeniz", x = names(all_dore_matrices))]
trojelsgaard_matrices <- all_dore_matrices[grepl("Trojelsgaard", x = names(all_dore_matrices))]
gibson_matrices <- all_dore_matrices[grepl("Gibson", x = names(all_dore_matrices))]

conditions <- expand.grid(
    model = c("iid", "pi", "rho", "pirho"),
    rep = seq.int(nb_rep),
    author = c("Baldock", "Traveset", "Souza", "Cordeniz", "Trojelsgaard", "Gibson")
)
conditions$model <- as.character(conditions$model)
conditions$author <- as.character(conditions$author)
set.seed(1234)
plan(list(
    tweak("callr", workers = floor(parallelly::availableCores(omit = 1L) / 3L)),
    tweak("callr", workers = 3L)
))
clustering_results <- future_lapply(seq_len(nrow(conditions)), function(s) {
    model <- conditions[s, "model"]
    rep <- conditions[s, "rep"]
    author <- conditions[s, "author"]

    if (author == "Baldock") {
        subdore_matrices <- baldock_matrices
    } else if (author == "Traveset") {
        subdore_matrices <- traveset_matrices
    } else if (author == "Souza") {
        subdore_matrices <- souza_matrices
    } else if (author == "Cordeniz") {
        subdore_matrices <- cordeniz_matrices
    } else if (author == "Trojelsgaard") {
        subdore_matrices <- trojelsgaard_matrices
    } else if (author == "Gibson") {
        subdore_matrices <- gibson_matrices
    }

    message(
        "Starting condition ", s, " on ", nrow(conditions),
        " with model = ", model, " author ", author, " and rep = ", rep
    )

    out <- clusterize_bipartite_networks(
        netlist = subdore_matrices,
        net_id = names(subdore_matrices),
        colsbm_model = model,
        global_opts = list(backend = "future"),
        fit_opts = list(max_vem_steps = 10000L)
    )

    message("End clustering ", s, " on ", nrow(conditions))
    saveRDS(out, file = file.path(temp_path, paste0("condition", s, "_on_", nrow(conditions), "_", model, "_", author, "_", rep, ".Rds")))
    return(out)
}, future.seed = TRUE)

saveRDS(clustering_results, file = file.path(save_path, paste0("clustering_results_", start_time, ".Rds")))
