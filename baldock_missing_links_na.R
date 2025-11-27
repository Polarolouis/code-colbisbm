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

# missing_links_list <- lapply(baldock_matrices, function(mat) {
#     nb_dyads_to_remove <- floor(max_prop_NAs * (nrow(mat) * ncol(mat)))

#     row_nodes <- sample.int(nrow(mat), size = nb_dyads_to_remove, replace = TRUE)
#     col_nodes <- sample.int(ncol(mat), size = nb_dyads_to_remove, replace = TRUE)

#     data.frame(
#         row = row_nodes,
#         col = col_nodes,
#         stringsAsFactors = FALSE
#     )
# })

missing_links_list <- lapply(baldock_matrices, function(mat) {
    links_df <- which(mat == 1, arr.ind = TRUE) |>
        as.data.frame() |>
        dplyr::rename(row = row, col = col)

    nb_links_to_remove <- floor(max_prop_NAs * nrow(links_df))
    sampled_indices <- sample.int(nrow(links_df), size = nb_links_to_remove, replace = FALSE)
    row_nodes <- links_df$row[sampled_indices]
    col_nodes <- links_df$col[sampled_indices]

    data.frame(
        row = row_nodes,
        col = col_nodes,
        stringsAsFactors = FALSE
    )
})

true_zeros_list <- lapply(baldock_matrices, function(mat) {
    zeros_df <- which(mat == 0, arr.ind = TRUE) |>
        as.data.frame() |>
        dplyr::rename(row = row, col = col)

    nb_zeros_to_remove <- floor(max_prop_NAs * nrow(zeros_df))
    sampled_indices <- sample.int(nrow(zeros_df), size = nb_zeros_to_remove, replace = FALSE)
    row_nodes <- zeros_df$row[sampled_indices]
    col_nodes <- zeros_df$col[sampled_indices]

    data.frame(
        row = row_nodes,
        col = col_nodes,
        stringsAsFactors = FALSE
    )
})

# For VGAE


epsilons <- seq(0.1, 0.8, by = 0.1)
possible_missing_network <- seq(1, length(baldock_matrices))
repetitions <- seq(1, 3)

vgae_conditions <- expand.grid(
    possible_missing_network = possible_missing_network,
    repetitions = repetitions,
    epsilon = epsilons
)
vgae_data_path <- file.path("data", "dore", "vgae_data")

write.csv(vgae_conditions, file.path(vgae_data_path, "baldock_missing_links_conditions.csv"), row.names = FALSE)


if (!dir.exists(vgae_data_path)) {
    dir.create(vgae_data_path)
}

lapply(
    seq_len(nrow(vgae_conditions)),
    function(s) {
        missing_network <- vgae_conditions[s, ]$possible_missing_network
        epsilon <- vgae_conditions[s, ]$epsilon
        repetition <- vgae_conditions[s, ]$repetitions

        missing_links_matrix <- baldock_matrices[[missing_network]]
        # names(missing_links_matrix) <- names(baldock_matrices)[missing_network]
        real_edge_label_index <- missing_links_list[[missing_network]]

        true_zeros_label_index <- true_zeros_list[[missing_network]]

        #  Remaining zeroes for train
        remaining_zeros <- true_zeros_label_index[!(true_zeros_label_index$row %in% real_edge_label_index$row & true_zeros_label_index$col %in% real_edge_label_index$col), ]


        # Selecting epsilon missing links
        real_edge_label_index <-
            real_edge_label_index[1:floor(epsilon * nrow(real_edge_label_index)), ]
        true_zeros_label_index <-
            true_zeros_label_index[1:nrow(real_edge_label_index), ]
        real_edge_label_index <- rbind(real_edge_label_index, true_zeros_label_index)
        real_values <- c()
        for (j in seq_len(nrow(real_edge_label_index))) {
            real_values <- c(real_values, missing_links_matrix[real_edge_label_index$row[j], real_edge_label_index$col[j]])
        }

        for (j in seq_len(nrow(real_edge_label_index))) {
            missing_links_matrix[real_edge_label_index$row[j], real_edge_label_index$col[j]] <- 0
        }

        train_edge_index <- which(missing_links_matrix == 1, arr.ind = TRUE)
        train_edge_index <- as.data.frame(train_edge_index)
        colnames(train_edge_index) <- c("row", "col")
        train_edge_index$split <- "train"
        train_edge_index$label <- 1

        train_zero_index <- remaining_zeros
        train_zero_index$split <- "train"
        train_zero_index$label <- 0


        real_edge_label_df <- real_edge_label_index
        real_edge_label_df$label <- real_values
        real_edge_label_df$split <- "test"

        combined_df <- rbind(train_edge_index, train_zero_index, real_edge_label_df)

        write.csv(combined_df, file.path(vgae_data_path, paste0("condition_", s, "_missing_network_", missing_network, "_epsilon_", epsilon, "_repetition_", repetition, ".csv")), row.names = FALSE)
    }
)

conditions <- expand.grid(
    possible_missing_network = possible_missing_network,
    repetitions = repetitions,
    epsilon = epsilons,
    model = c("iid", "pi", "rho", "pirho", "sep"),
    missing_replacement = c(NA, 0)
)

results <- future_lapply(seq_len(nrow(conditions)), function(s) {
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
    true_zeroes <- true_zeros_list[[missing_network]]
    true_zeroes <- true_zeroes[seq_len(nrow(missing_links)), ]
    missing_links <- rbind(missing_links, true_zeroes)
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
        netlist = ifelse(model == "sep", list(missing_links_matrix), matrices), # Only giving one matrix if we test sep model
        colsbm_model = model,
        net_id = names(matrices),
        nb_run = 1L,
        global_opts = list(backend = "no_mc")
    )

    tau1 <- fit$best_fit$tau[[missing_network]][[1]]
    alpha <- fit$best_fit$alpha
    tau2 <- fit$best_fit$tau[[missing_network]][[2]]

    predicted_values <- c()
    Xhat <- tau1 %*% alpha %*% t(tau2)
    for (j in seq_len(nrow(missing_links))) {
        predicted_values <- c(predicted_values, Xhat[missing_links$row[j], missing_links$col[j]])
    }

    # Compute ROC AUC for missings
    auc_value <- auc(real_values, predicted_values)
    out_df <- cbind(conditions[s, ], data.frame(AUC = auc_value))
    saveRDS(out_df, file = file.path(temp_path, paste0("condition_", s, "_on_", nrow(conditions), ".Rds")))
    return(out_df)
}, future.seed = TRUE)

saveRDS(do.call("rbind", results), file = save_file)
