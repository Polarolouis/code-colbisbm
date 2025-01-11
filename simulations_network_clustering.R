necessary_packages <- c("remotes", "colSBM")
library(here)

options(future.globals.maxSize = Inf)

if (!all(necessary_packages %in% installed.packages())) {
    install.packages(necessary_packages[-length(necessary_packages)])
}
suppressPackageStartupMessages(library("colSBM"))
future::plan(future::multisession())

set.seed(0L)


nr <- 75
nc <- 75

pi <- matrix(c(0.2, 0.3, 0.5), nrow = 1, byrow = TRUE)
rho <- matrix(c(0.2, 0.3, 0.5), nrow = 1, byrow = TRUE)
repetitions <- seq.int(30)
epsilons <- seq(0.1, 0.4, by = 0.1)
models <- c("iid", "pi", "rho", "pirho")

save_folder <- here(
    "simulations", "clustering",
    "9collection"
)

if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
}

save_filename <- paste0(
    "9collection_data_clustering_",
    format(Sys.time(), "%d-%m-%y-%H-%M-%S"),
    ".Rds"
)

temp_folder <- file.path(save_folder, paste0("tmp", format(Sys.time(), "%d-%m-%y-%H-%M-%S")))

if (!dir.exists(temp_folder)) {
    dir.create(temp_folder, recursive = TRUE)
}

conditions <- tidyr::crossing(epsilons, repetitions, models)
with_progress({
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

            list_collection <- clusterize_bipartite_networks_graphon(
                netlist = incidence_matrices,
                net_id = netids,
                nb_run = 3L,
                colsbm_model = current_model,
                global_opts = list(
                    nb_cores = parallelly::availableCores(omit = 1L), verbosity = 2,
                    plot_details = 0 # ,
                    # parallelization_vector = c(FALSE, FALSE, FALSE)
                )
            )

            best_partitions <- unlist(tail(list_collection, 1))
            if (!is(best_partitions, "list")) {
                best_partitions <- list(best_partitions)
            }
            clustering <- unlist(lapply(seq_along(best_partitions), function(col_idx) {
                setNames(
                    rep(col_idx, best_partitions[[col_idx]]$M),
                    best_partitions[[col_idx]]$net_id
                )
            }))
            # ARI computation
            clustering <- clustering[order(names(clustering))]
            ari <- aricode::ARI(rep(c(1, 2, 3), each = 3), clustering)
            out <- data.frame(epsilon = eps, model = current_model, ARI = ari)

            saveRDS(out, file = file.path(
                temp_folder,
                paste0("condition_", s, "_on_", nrow(conditions), ".Rds")
            ))

            return(out)
        },
        future.seed = NULL
    )
})
data_frame_result <- do.call("rbind", results)

saveRDS(data_frame_result, file = file.path(
    save_folder,
    save_filename
))
