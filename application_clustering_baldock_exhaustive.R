library(colSBM)
library(ConsRank)
library(future)
library(future.apply)
library(future.callr)
library(here)

options(future.globals.maxSize = Inf)

(all_partitions <- partitions(5))
(all_partitions <- all_partitions[seq(nrow(all_partitions), 1), ])
# Here we want to merge the strings for each row e.g "{a", "b"}" into "ab"

partitions_keys <- lapply(seq_len(nrow(all_partitions)), function(idx) {
    # Here we merge the strings if there is an opening { and until a closing }
    # is met we also remove the opening and closing brackets

    curr_string <- all_partitions[idx, ]
    new_string <- curr_string
    merge <- FALSE
    merge_idx <- 1
    for (i in seq_along(curr_string)) {
        if (merge) {
            new_string[merge_idx] <- paste0(new_string[merge_idx], curr_string[i])
            new_string[i] <- ""
            if (grepl("\\}", curr_string[i])) {
                merge <- FALSE
                new_string[merge_idx] <- gsub("\\}", "", new_string[merge_idx])
            }
        } else {
            new_string[i] <- curr_string[i]
            if (grepl("\\{", curr_string[i])) {
                merge <- TRUE
                merge_idx <- i
                new_string[i] <- gsub("\\{", "", curr_string[i])
            }
        }
    }

    new_string <- new_string[new_string != ""]
})

baldock_matrices <- readRDS(here("data", "dore-binary-matrices.Rds"))
baldock_matrices <- baldock_matrices[grepl("Baldock", names(baldock_matrices))]

match_string_to_netids <- function(string, netids = names(baldock_matrices)) {
    letters_to_match <- unlist(strsplit(string, ""))

    sapply(letters_to_match, function(letter) {
        switch(letter,
            "a" = netids[1],
            "b" = netids[2],
            "c" = netids[3],
            "d" = netids[4],
            "e" = netids[5]
        )
    }) |> unname()
}

memoization_list <- list()

is_memoized <- function(key) {
    # Check if the key is already in memoization_list
    if (key %in% names(memoization_list)) {
        return(TRUE)
    }
    return(FALSE)
}

compute_collection <- function(key) {
    netids <- match_string_to_netids(key)
    matrices <- baldock_matrices[netids]
    fit <- estimate_colBiSBM(
        netlist = matrices,
        net_id = netids,
        colsbm_model = "iid",
        global_opts = list(backend = "future"),
        fit_opts = list(max_vem_steps = 10000L)
    )

    return(fit)
}

compute_partition <- function(partition_keys, netids = names(baldock_matrices)) {
    memoized_bool <- sapply(partition_keys, is_memoized)
    not_memoized_keys <- partition_keys[!memoized_bool]

    fits_not_memoized <- future_lapply(not_memoized_keys, function(key) {
        fit <- compute_collection(key)
        return(fit)
    }, future.seed = TRUE)
    names(fits_not_memoized) <- not_memoized_keys
    # Store the results in the memoization list
    memoization_list <<- c(memoization_list, fits_not_memoized)

    #  Memoized results are already in the list
    results <- memoization_list[partition_keys[memoized_bool]]
    results <- c(results, fits_not_memoized)
    results <- results[partition_keys]
    return(results)
}

plan(list(tweak("callr", workers = 5L), tweak("callr", workers = 3L)))

save_path <- here("applications", "baldock_exhaustive")

if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
}

start_time <- format(Sys.time(), "%Y%m%d%H%M%S")
exhaustive_partitions_list <- lapply(seq_len(nrow(all_partitions)), function(idx) {
    message("Starting partition ", idx, " of ", nrow(all_partitions))
    partition_keys <- partitions_keys[[idx]]
    out <- compute_partition(partition_keys)
    saveRDS(out, file = file.path(save_path, paste0(start_time, "exhaustive_partition_", idx, "_on_", nrow(all_partitions), ".Rds")))
})
message("Finished all partitions")
