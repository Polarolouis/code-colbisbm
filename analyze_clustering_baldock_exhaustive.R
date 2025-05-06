library(colSBM)
library(here)
library(stringr)
library(ggplot2)
library(patchwork)
devtools::load_all("../colSBM/")

iid_data_files <- list.files(here("applications", "baldock_exhaustive"), full.names = TRUE) |> str_sort(numeric = TRUE)

extract_best_partition_from_exh <- function(data_files) {
    all_partitions <- lapply(data_files, readRDS)
    all_bicls <- sapply(all_partitions, compute_bicl_partition)
    all_partitions <- all_partitions[order(all_bicls, decreasing = TRUE)]
    all_partitions <- lapply(all_partitions, function(partition) {
        lapply(partition, function(col) {
            col$best_fit
        })
    })
    all_bicls <- all_bicls[order(all_bicls, decreasing = TRUE)]
    best_partition <- all_partitions[[which.max(all_bicls)]]

    # lapply(all_partitions, function(best_partition) {
    #     print(wrap_plots(
    #         {
    #             lapply(best_partition, function(col) {
    #                 free(plot(col, type = "meso", values = TRUE, mixture = TRUE), type = "label")
    #             })
    #         },
    #         nrow = 2
    #     ) +
    #         plot_layout(guides = "collect") +
    #         plot_annotation(title = "Best partition") & theme(plot.title = element_text(hjust = 0.5), legend.position = "right"))
    # })
    return(best_partition)
}

pirho_data_files <- list.files(here("applications", "baldock_exhaustive", "pirho"), full.names = TRUE) |> str_sort(numeric = TRUE)

# extract_best_partition_from_exh(iid_data_files)
