necessary_packages <- c("remotes", "colSBM", "future.apply", "future.callr", "progressr")

options(future.globals.maxSize = Inf)

if (!all(necessary_packages %in% installed.packages())) {
    install.packages(necessary_packages)
}
library(here)
library(future.apply)
library(future.callr)
library(progressr)
suppressPackageStartupMessages(library("colSBM"))
handlers(global = TRUE)
plan(list(
    tweak(callr, workers = parallelly::availableCores(omit = 2L) / 9L),
    tweak(callr, workers = 9L)
))

set.seed(0L)


nr <- 75
nc <- 75

pi <- matrix(c(0.05, 0.3, 0.65), nrow = 1, byrow = TRUE)
rho <- matrix(c(0.1, 0.8, 0.1), nrow = 1, byrow = TRUE)
repetitions <- seq.int(10)
epsilons <- seq(0.1, 0.4, by = 0.1)
models <- c("iid") # , "pi", "rho", "pirho" not implemented

save_folder <- here(
    "simulations", "clustering",
    "9collection"
)

if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
}

save_filename <- paste0(
    "9collection_data_clustering_asc",
    format(Sys.time(), "%d-%m-%y-%H-%M-%S"),
    ".Rds"
)

temp_folder <- file.path(save_folder, paste0("tmp", format(Sys.time(), "%d-%m-%y-%H-%M-%S")))

if (!dir.exists(temp_folder)) {
    dir.create(temp_folder, recursive = TRUE)
}

conditions <- tidyr::crossing(epsilons, repetitions, models)
with_progress(
    {
        pb <- progressr::progressor(along = seq_len(nrow(conditions)))

        results <- future.apply::future_lapply(
            seq_len(nrow(conditions)), function(s) {
                eps <- conditions[s, ]$epsilons
                current_pi <- pi
                current_rho <- rho
                current_model <- conditions[s, ]$models

                cli::cli_text("Starting condition {s} on {nrow(conditions)}")

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
                    current_pi, current_rho,
                    alpha_assortative, 3,
                    model = current_model,
                    return_memberships = TRUE
                )

                assortative_incidence <- lapply(
                    seq_along(assortative_collection),
                    function(m) {
                        return(assortative_collection[[m]]$incidence_matrix)
                    }
                )

                assortative_row_clustering <- lapply(
                    seq_along(assortative_collection),
                    function(m) {
                        return(assortative_collection[[m]]$row_clustering)
                    }
                )

                assortative_col_clustering <- lapply(
                    seq_along(assortative_collection),
                    function(m) {
                        return(assortative_collection[[m]]$row_clustering)
                    }
                )

                core_periphery_collection <- generate_bipartite_collection(
                    nr, nc,
                    current_pi, current_rho,
                    alpha_core_periphery, 3,
                    model = current_model,
                    return_memberships = TRUE
                )

                core_periphery_incidence <- lapply(
                    seq_along(core_periphery_collection),
                    function(m) {
                        return(core_periphery_collection[[m]]$incidence_matrix)
                    }
                )

                core_periphery_row_clustering <- lapply(
                    seq_along(core_periphery_collection),
                    function(m) {
                        return(core_periphery_collection[[m]]$row_clustering)
                    }
                )

                core_periphery_col_clustering <- lapply(
                    seq_along(core_periphery_collection),
                    function(m) {
                        return(core_periphery_collection[[m]]$row_clustering)
                    }
                )

                disassortative_collection <- generate_bipartite_collection(
                    nr, nc,
                    current_pi, current_rho,
                    alpha_disassortative, 3,
                    model = current_model,
                    return_memberships = TRUE
                )

                disassortative_incidence <- lapply(
                    seq_along(disassortative_collection),
                    function(m) {
                        return(disassortative_collection[[m]]$incidence_matrix)
                    }
                )

                disassortative_row_clustering <- lapply(
                    seq_along(disassortative_collection),
                    function(m) {
                        return(disassortative_collection[[m]]$row_clustering)
                    }
                )

                disassortative_col_clustering <- lapply(
                    seq_along(disassortative_collection),
                    function(m) {
                        return(disassortative_collection[[m]]$row_clustering)
                    }
                )

                real_row_clustering <- append(
                    append(
                        assortative_row_clustering,
                        core_periphery_row_clustering
                    ),
                    disassortative_row_clustering
                )

                real_col_clustering <- append(
                    append(
                        assortative_col_clustering,
                        core_periphery_col_clustering
                    ),
                    disassortative_col_clustering
                )

                incidence_matrices <- append(
                    append(
                        assortative_incidence,
                        core_periphery_incidence
                    ),
                    disassortative_incidence
                )

                netids <- rep(c("as", "cp", "dis"), each = 3)

                clust_out <- clusterize_bipartite_networks_graphon(
                    netlist = incidence_matrices,
                    net_id = netids,
                    nb_run = 3L,
                    colsbm_model = current_model,
                    global_opts = list(
                        backend = "future",
                        verbosity = 1L
                    ),
                    fusions_per_step = 1L,
                    fit_opts = list(max_vem_steps = 3000L)
                )

                best_partitions <- tail(clust_out$fusion_history, 1)[[1]]
                bicl_vec <- sapply(seq_along(best_partitions), function(col_idx) {
                    return(best_partitions[[col_idx]]$best_fit$BICL)
                })
                bicl <- ifelse(length(bicl_vec) > 1, sum(bicl_vec), ifelse(length(bicl_vec) == 1, bicl_vec, NA))
                clustering <- unlist(lapply(seq_along(best_partitions), function(col_idx) {
                    setNames(
                        rep(col_idx, best_partitions[[col_idx]]$M),
                        best_partitions[[col_idx]]$net_id
                    )
                }))
                # ARI computation
                clustering <- clustering[order(names(clustering))]
                ari <- aricode::ARI(rep(c(1, 2, 3), each = 3), clustering)
                out <- data.frame(epsilon = eps, model = current_model, ARI = ari, BICL = bicl)

                saveRDS(out, file = file.path(
                    temp_folder,
                    paste0("condition_", s, "_on_", nrow(conditions), ".Rds")
                ))

                return(out)
            },
            future.seed = NULL
        )
    },
    delay_stdout = TRUE,
    delay_conditions = "condition"
)
data_frame_result <- do.call("rbind", results)

saveRDS(data_frame_result, file = file.path(
    save_folder,
    save_filename
))
