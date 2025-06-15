necessary_packages <- c("remotes", "colSBM", "future.apply", "future.callr", "progressr")

options(future.globals.maxSize = Inf)

stopifnot("Function `clusterize_bipartite_networks_d_a` not found" = tryCatch(
    {
        getFromNamespace(x = "clusterize_bipartite_networks_d_a", ns = "colSBM")
        TRUE
    },
    error = function(...) FALSE
))

if (!all(necessary_packages %in% installed.packages())) {
    install.packages(necessary_packages)
}
library(here)
library(future)
library(future.apply)
# library(future.callr)
library(progressr)
suppressPackageStartupMessages(library("colSBM"))
handlers(global = TRUE)
plan(list(
    tweak("multisession", workers = parallelly::availableCores(omit = 2L) / 3L),
    tweak("multisession", workers = 3L)
))

set.seed(0L)


nr <- 75
nc <- 75
M <- 3L

pi <- matrix(c(0.05, 0.3, 0.65), nrow = 1, byrow = TRUE)
rho <- matrix(c(0.1, 0.8, 0.1), nrow = 1, byrow = TRUE)
repetitions <- seq.int(10)
epsilons <- seq(0.1, 0.4, by = 0.1)
models <- c("iid")

save_folder <- here(
    "simulations", "clustering",
    paste0(3 * M, "collection")
)

if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
}

save_filename <- paste0(
    paste0(3 * M, "collection_data_clustering_da_"),
    format(Sys.time(), "%d-%m-%y-%H-%M-%S"),
    ".Rds"
)

temp_folder <- file.path(save_folder, paste0("tmp", format(Sys.time(), "%d-%m-%y-%H-%M-%S")))

if (!dir.exists(temp_folder)) {
    dir.create(temp_folder, recursive = TRUE)
}

conditions <- tidyr::crossing(epsilons, repetitions, models)

# missing_conditions_file <- file.path(save_folder, "missing_conditions_da.Rds")

# if (file.exists(missing_conditions_file)) {
#     message("Resuming from missing conditions.")
#     row_conditions <- readRDS(missing_conditions_file)
# } else {
message("Starting from scratch.")
row_conditions <- seq_len(nrow(conditions))
# }

# Pre-generate the matrices list to have good comparison
matrices_lists <- lapply(models, function(model) {
    lapply(epsilons, function(eps) {
        lapply(repetitions, function(rep) {
            alpha_assortative <- matrix(0.3, nrow = 3, ncol = 3) +
                matrix(
                    c(
                        eps, -0.5 * eps, -0.5 * eps,
                        -0.5 * eps, eps, -0.5 * eps,
                        -0.5 * eps, -0.5 * eps, eps
                    ),
                    nrow = 3, byrow = TRUE
                )

            alpha_core_periphery <- matrix(0.3, nrow = 3, ncol = 3) +
                matrix(
                    c(
                        1.5 * eps, eps, 0.5 * eps,
                        eps, 0.5 * eps, 0,
                        0.5 * eps, 0, -0.5 * eps
                    ),
                    nrow = 3, byrow = TRUE
                )

            alpha_disassortative <- matrix(0.3, nrow = 3, ncol = 3) +
                matrix(
                    c(
                        -0.5 * eps, eps, eps,
                        eps, -0.5 * eps, eps,
                        eps, eps, -0.5 * eps
                    ),
                    nrow = 3, byrow = TRUE
                )

            assortative_collection <- generate_bipartite_collection(
                nr, nc,
                pi, rho,
                alpha_assortative, M,
                model = model,
                return_memberships = TRUE
            )

            assortative_incidence <- lapply(
                seq_along(assortative_collection),
                function(m) {
                    return(assortative_collection[[m]]$incidence_matrix)
                }
            )

            core_periphery_collection <- generate_bipartite_collection(
                nr, nc,
                pi, rho,
                alpha_core_periphery, M,
                model = model,
                return_memberships = TRUE
            )

            core_periphery_incidence <- lapply(
                seq_along(core_periphery_collection),
                function(m) {
                    return(core_periphery_collection[[m]]$incidence_matrix)
                }
            )

            disassortative_collection <- generate_bipartite_collection(
                nr, nc,
                pi, rho,
                alpha_disassortative, M,
                model = model,
                return_memberships = TRUE
            )

            disassortative_incidence <- lapply(
                seq_along(disassortative_collection),
                function(m) {
                    return(disassortative_collection[[m]]$incidence_matrix)
                }
            )


            incidence_matrices <- append(
                append(
                    assortative_incidence,
                    core_periphery_incidence
                ),
                disassortative_incidence
            )
            netids <- paste0(rep(c("as", "cp", "dis"), each = M), ".", seq(1, M))
            incidence_matrices <- setNames(incidence_matrices, netids)

            return(incidence_matrices)
        })
    })
})

names(matrices_lists) <- models
for (model in models) {
    names(matrices_lists[[model]]) <- epsilons
    for (eps in epsilons) {
        names(matrices_lists[[model]][[as.character(eps)]]) <- repetitions
    }
}

results <- future.apply::future_lapply(
    seq_len(nrow(conditions)), function(s) {
        if (!(s %in% row_conditions)) {
            message("Skipping condition ", s, " on ", nrow(conditions))
            return(NULL)
        }
        eps <- conditions[s, ]$epsilons
        current_pi <- pi
        current_rho <- rho
        current_model <- conditions[s, ]$models
        rep <- conditions[s, ]$repetitions

        cli::cli_text("Starting condition {s} on {nrow(conditions)}")


        netids <- paste0(rep(c("as", "cp", "dis"), each = M), ".", seq(1, M))
        incidence_matrices <- matrices_lists[[current_model]][[as.character(eps)]][[rep]]

        list_collection <- clusterize_bipartite_networks_d_a(
            netlist = incidence_matrices,
            net_id = netids,
            nb_run = 3L,
            colsbm_model = current_model,
            global_opts = list(
                backend = "future",
                verbosity = 1L
            ),
            fit_opts = list(max_vem_steps = 10000L)
        )

        best_partitions <- list_collection$partition
        bicl_vec <- sapply(seq_along(best_partitions), function(col_idx) {
            return(best_partitions[[col_idx]]$BICL)
        })
        bicl <- ifelse(length(bicl_vec) > 1, sum(bicl_vec), ifelse(length(bicl_vec) == 1, bicl_vec, NA))
        # ARI computation
        clustering <- list_collection$cluster
        ari <- aricode::ARI(rep(c(1, 2, 3), each = M), clustering)
        out <- data.frame(epsilon = eps, model = current_model, ARI = ari, BICL = bicl, nb_collections = length(best_partitions))

        saveRDS(out, file = file.path(
            temp_folder,
            paste0("condition_", s, "_on_", nrow(conditions), ".Rds")
        ))

        return(out)
    },
    future.seed = TRUE
)

data_frame_result <- do.call("rbind", results)

saveRDS(data_frame_result, file = file.path(
    save_folder,
    save_filename
))
