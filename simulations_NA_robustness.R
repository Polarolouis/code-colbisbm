necessary_packages <- c(
    "remotes", "tictoc", "pROC", "parallel", "colSBM", "future.apply"
)

options(future.globals.maxSize = Inf)

library("pROC")
library("colSBM")

set.seed(1234)

eps <- 0.05

M <- 2

# Defining parameters
nr1 <- 20L
nc1 <- 20L
nr2 <- 120
nc2 <- 120
pir <- c(0.5, 0.3, 0.2)
pic <- c(0.5, 0.3, 0.2)

alpha_modular <- matrix(c(
    0.9, eps, eps,
    eps, 0.2, eps,
    eps, eps, 0.8
), byrow = TRUE, nrow = length(pir), ncol = length(pic))

alpha_nested <- matrix(c(
    0.9, 0.65, 0.1,
    0.35, 0.15, eps,
    0.1, eps, eps
), byrow = TRUE, nrow = length(pir), ncol = length(pic))

max_repetition <- 10L

# Collections
collections <- list(
    "iid" = list("modular" = c(
        generate_bipartite_collection(nr1, nc1,
            pir, pic,
            alpha_modular,
            M = 1,
            model = "iid",
            return_memberships = TRUE
        ),
        generate_bipartite_collection(nr2, nc2,
            pir, pic,
            alpha_modular,
            M = M - 1,
            model = "iid",
            return_memberships = TRUE
        )
    ), "nested" = c(
        generate_bipartite_collection(nr1, nc1,
            pir, pic,
            alpha_nested,
            M = 1,
            model = "iid",
            return_memberships = TRUE
        ),
        generate_bipartite_collection(nr2, nc2,
            pir, pic,
            alpha_nested,
            M = M - 1,
            model = "iid",
            return_memberships = TRUE
        )
    )),
    "pi" = list(
        "modular" = c(
            generate_bipartite_collection(nr1, nc1,
                pir, pic,
                alpha_modular, 1,
                model = "pi",
                return_memberships = TRUE
            ),
            generate_bipartite_collection(nr2, nc2,
                pir, pic,
                alpha_modular, M - 1,
                model = "pi",
                return_memberships = TRUE
            )
        ),
        "nested" = c(
            generate_bipartite_collection(nr1, nc1,
                pir, pic,
                alpha_nested, 1,
                model = "pi",
                return_memberships = TRUE
            ),
            generate_bipartite_collection(nr2, nc2,
                pir, pic,
                alpha_nested, M - 1,
                model = "pi",
                return_memberships = TRUE
            )
        )
    ),
    "rho" = list(
        "modular" = c(
            generate_bipartite_collection(nr1, nc1,
                pir, pic,
                alpha_modular, 1,
                model = "rho",
                return_memberships = TRUE
            ),
            generate_bipartite_collection(nr2, nc2,
                pir, pic,
                alpha_modular, M - 1,
                model = "rho",
                return_memberships = TRUE
            )
        ),
        "nested" = c(
            generate_bipartite_collection(nr1, nc1,
                pir, pic,
                alpha_nested, 1,
                model = "rho",
                return_memberships = TRUE
            ),
            generate_bipartite_collection(nr2, nc2,
                pir, pic,
                alpha_nested, M - 1,
                model = "rho",
                return_memberships = TRUE
            )
        )
    ),
    "pirho" = list(
        "modular" = c(
            generate_bipartite_collection(nr1, nc1,
                pir, pic,
                alpha_modular, 1,
                model = "pirho",
                return_memberships = TRUE
            ),
            generate_bipartite_collection(nr2, nc2,
                pir, pic,
                alpha_modular, M - 1,
                model = "pirho",
                return_memberships = TRUE
            )
        ),
        "nested" = c(
            generate_bipartite_collection(nr1, nc1,
                pir, pic,
                alpha_nested, 1,
                model = "pirho",
                return_memberships = TRUE
            ),
            generate_bipartite_collection(nr2, nc2,
                pir, pic,
                alpha_nested, M - 1,
                model = "pirho",
                return_memberships = TRUE
            )
        )
    )
)


conditions <- expand.grid(
    prop_NAs = seq(from = 0, to = 0.9, by = 0.1),
    model = c("iid", "pi", "rho", "pirho"),
    struct = c("modular", "nested"),
    repetition = seq.int(max_repetition)
)


#  Data params
main_dir <- file.path("simulations", "NA_robustness")

if (!dir.exists(main_dir)) {
    dir.create(main_dir, recursive = TRUE)
}

start_time <- format(Sys.time(), "%d-%m-%Y_%H-%M-%S")
temp_dir <- file.path(main_dir, paste0(
    "tmp", start_time
))

if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
}

file_save <- file.path(main_dir, paste0(
    "NA_robustness_", start_time, "_1-", nrow(conditions), ".Rds"
))

message(
    "Starting NA robustness simulation."
)

library("future")
library("future.callr")
plan(
    list(
        tweak("callr", workers = floor(parallelly::availableCores(omit = 1L) / 3L)),
        tweak("callr", workers = 3L)
    )
)
result_list <- future.apply::future_lapply(
    seq_len(nrow(conditions)),
    function(current) {
        # Looping over conditions
        prop_NAs <- conditions[current, ][["prop_NAs"]]
        model <- as.character(conditions[current, ][["model"]])
        struct <- as.character(conditions[current, ][["struct"]])
        bipartite_collection <- collections[[model]][[struct]]

        # This is a list of the M incidence matrices
        bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
            bipartite_collection[[m]][["incidence_matrix"]]
        })

        Z <- lapply(seq.int(M), function(m) {
            list(
                bipartite_collection[[m]][["row_blockmemberships"]],
                bipartite_collection[[m]][["col_blockmemberships"]]
            )
        })

        NAs_selected_index <- seq_len(length(bipartite_collection_incidence[[1]]))

        NAs_index <- sample(NAs_selected_index, size = floor(prop_NAs * length(NAs_selected_index)))

        real_val_NAs <- bipartite_collection_incidence[[1]][NAs_index]
        bipartite_collection_incidence[[1]][NAs_index] <- NA
        NAs_coordinates <- which(is.na(bipartite_collection_incidence[[1]]),
            arr.ind = TRUE
        )

        start_time <- Sys.time()
        mybisbmpop <- estimate_colBiSBM(
            netlist = bipartite_collection_incidence, colsbm_model = model,
            nb_run = 3L,
            global_opts = list(
                nb_cores = parallel::detectCores() - 1, verbosity = 0,
                backend = "future"
            )
        )
        stop_time <- Sys.time()

        baseline_LBM <- estimate_colBiSBM(
            netlist = list(bipartite_collection_incidence[[1]]), colsbm_model = "iid",
            nb_run = 3L,
            global_opts = list(
                nb_cores = parallel::detectCores() - 1, verbosity = 0,
                backend = "future"
            )
        )

        # Predicted links
        X_hat_LBM <- baseline_LBM[["best_fit"]][["tau"]][[1]][[1]] %*%
            baseline_LBM[["best_fit"]][["alpha"]] %*%
            t(baseline_LBM[["best_fit"]][["tau"]][[1]][[2]])
        X_hat <- mybisbmpop[["best_fit"]][["tau"]][[1]][[1]] %*%
            mybisbmpop[["best_fit"]][["alpha"]] %*%
            t(mybisbmpop[["best_fit"]][["tau"]][[1]][[2]])

        # Compute ROC and AUC
        auc_LBM <- auc(c(0, 1, real_val_NAs), c(0, 1, X_hat_LBM[NAs_index]))
        auc_colBiSBM <- auc(c(0, 1, real_val_NAs), c(0, 1, X_hat[NAs_index]))

        # Computing ARI on the NAs
        out_data_frame <- data.frame(
            prop_NAs = prop_NAs,
            model = model,
            struct = struct,
            repetition = conditions[current, ][["repetition"]],
            auc_LBM = auc_LBM,
            auc_colBiSBM = auc_colBiSBM,
            arirow_LBM = aricode::ARI(
                Z[[1]][[1]],
                baseline_LBM[["best_fit"]][["Z"]][[1]][[1]]
            ),
            aricol_LBM = aricode::ARI(
                Z[[1]][[2]],
                baseline_LBM[["best_fit"]][["Z"]][[1]][[2]]
            ),
            arirow_colBiSBM = aricode::ARI(
                Z[[1]][[1]],
                mybisbmpop[["best_fit"]][["Z"]][[1]][[1]]
            ),
            aricol_colBiSBM = aricode::ARI(
                Z[[1]][[2]],
                mybisbmpop[["best_fit"]][["Z"]][[1]][[2]]
            ),
            elapsed_secs = difftime(stop_time, start_time, units = "sec")
        )

        message("Finished step ", current, "/", nrow(conditions), "\n")

        #  Saving temp
        temp_file_save <- file.path(temp_dir, paste0(
            "conditions_", current, "_on_",
            nrow(conditions), ".Rds"
        ))

        saveRDS(out_data_frame, file = temp_file_save)

        return(out_data_frame)
    },
    future.seed = NULL
)

result_dataframe <- do.call("rbind", result_list)


saveRDS(
    result_dataframe,
    file = file_save
)
message("Finished simulations.")
