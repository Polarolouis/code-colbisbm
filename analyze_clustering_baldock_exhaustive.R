library(colSBM)
library(here)
library(stringr)
library(ggplot2)
library(patchwork)

data_files <- list.files(here("applications", "baldock_exhaustive"), full.names = TRUE) |> str_sort(numeric = TRUE)

all_partitions <- lapply(data_files, readRDS)
all_bicls <- sapply(all_partitions, compute_bicl_partition)

best_partition <- all_partitions[[which.max(all_bicls)]]

wrap_plots(
    {
        lapply(best_partition, function(col) {
            free(plot(col$best_fit, type = "meso", values = TRUE, mixture = TRUE), type = "label")
        })
    },
    nrow = 2
) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Best partition") & theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
