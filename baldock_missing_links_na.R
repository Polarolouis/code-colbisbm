library(colSBM)
library(here)
library(pROC)
library(future.apply)

plan("multisession")
epoch <- as.numeric(as.POSIXct(Sys.time(), tz = "GMT", origin = "1970-01-01"))
save_path <- here("simulations", "missing_links")
save_file <- file.path(save_path, paste0("missing_links_baldock_", epoch, ".Rds"))
temp_path <- file.path(save_path, paste0("tmp", epoch))
if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
}
if (!dir.exists(temp_path)) {
    dir.create(temp_path, recursive = TRUE)
}
set.seed(42)
baldock_matrices <- readRDS(here("data", "baldock2019-binary-matrices.Rds"))

max_prop_NAs <- 0.9

missing_links_list <- lapply(baldock_matrices, function(mat) {
    nb_links_to_remove <- floor(max_prop_NAs * (nrow(mat) * ncol(mat)))

    row_nodes <- sample.int(nrow(mat), size = nb_links_to_remove, replace = TRUE)
    col_nodes <- sample.int(ncol(mat), size = nb_links_to_remove, replace = TRUE)

    data.frame(
        row = row_nodes,
        col = col_nodes,
        stringsAsFactors = FALSE
    )
})

epsilons <- seq(0.1, 0.8, by = 0.1)
possible_missing_network <- seq(1, length(baldock_matrices))
repetitions <- seq(1, 3)

conditions <- expand.grid(
    possible_missing_network = possible_missing_network,
    repetitions = repetitions,
    epsilon = epsilons,
    model = c("iid", "pi", "rho", "pirho"),
    missing_replacement = c(NA, 0)
)

results <- future_lapply(seq_along(conditions), function(s) {
    missing_network <- conditions[s, ]$possible_missing_network
    epsilon <- conditions[s, ]$epsilon
    repetition <- conditions[s, ]$repetitions
    model <- as.character(conditions[s, ]$model)
    missing_replacement <- conditions[s, ]$missing_replacement

    message(
        sprintf(
            "Missing network %s, epsilon %.1f, repetition %d, model %s, missing replacement %s",
            names(baldock_matrices)[missing_network], epsilon, repetition, model,
            ifelse(is.na(missing_replacement), "NA", missing_replacement)
        )
    )

    complete_matrices <- baldock_matrices[-missing_network]
    missing_links_matrix <- baldock_matrices[[missing_network]]
    names(missing_links_matrix) <- names(baldock_matrices)[missing_network]
    missing_links <- missing_links_list[[missing_network]]
    # Selecting epsilon missing links
    missing_links <- missing_links[1:floor(epsilon * nrow(missing_links)), ]
    real_values <- c()
    for (j in seq_len(nrow(missing_links))) {
        real_values <- c(real_values, missing_links_matrix[missing_links$row[j], missing_links$col[j]])
        missing_links_matrix[missing_links$row[j], missing_links$col[j]] <- missing_replacement
    }

    matrices <- append(
        complete_matrices, setNames(list(missing_links_matrix), names(baldock_matrices)[missing_network]),
        after = missing_network - 1
    )

    fit <- colSBM::estimate_colBiSBM(
        netlist = matrices,
        colsbm_model = model,
        net_id = names(matrices),
        nb_run = 1L,
        global_opts = list(backend = "no_mc")
    )

    predicted_values <- c()
    Xhat <- fit$best_fit$tau[[missing_network]][[1]] %*% fit$best_fit$alpha %*% t(fit$best_fit$tau[[missing_network]][[2]])
    for (j in seq_len(nrow(missing_links))) {
        predicted_values <- c(predicted_values, Xhat[missing_links$row[j], missing_links$col[j]])
    }

    # Compute ROC AUC for missings
    auc_value <- auc(real_values, predicted_values)
    out_df <- cbind(conditions[s, ], data.frame(AUC = auc_value))
    saveRDS(out_df, file = file.path(temp_path, paste0("condition_", s, "_on_", nrow(conditions), ".Rds")))
    return(out_df)
})
