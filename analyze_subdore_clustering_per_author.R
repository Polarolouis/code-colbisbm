library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ggalluvial)
library(ggokabeito)

data_path <- here("applications", "subdore", "tmp20250418103624")

author_vector <- c("Baldock", "Traveset", "Souza", "Cordeniz", "Trojelsgaard", "Gibson")
model_vector <- c("iid", "pi", "rho", "pirho")


flist <- list.files(data_path)
author_order <- order(factor(sub(".*_(Baldock|Traveset|Souza|Cordeniz|Trojelsgaard|Gibson)_.*", "\\1", flist), levels = author_vector))
flist[author_order]

# Create named nested list
# $Author$iid
# $Author$pi
# $Author$rho
# $Author$pirho

# setNames(lapply(author_vector, function(author) {
#     # Filter flist by author
#     author_files <- flist[grepl(author, flist)]
#     # Filte by model
#     setNames(lapply(model_vector, function(model) {
#         # Filter author_files by model
#         model_files <- author_files[grepl(model, author_files)]
#         # Read files
#         lapply(model_files, function(file) {
#             readRDS(file.path(data_path, file))
#         })
#     }), paste0(model_vector, seq(1, 4)))
# }), author_vector) -> nested_list

baldock_flist <- flist[grepl("Baldock", flist)]
baldock_flist <- baldock_flist[order(factor(sub(".*_([1-5]).Rds", "\\1", baldock_flist), levels = seq(1, 5)))]
baldock_flist <- baldock_flist[order(factor(sub(".*_(iid|pi|rho|pirho)_.*", "\\1", baldock_flist), levels = model_vector))]
baldock_flist

extract_from_flist <- function(flist) {
    setNames(lapply(model_vector, function(model) {
        # Filter author_files by model
        model_files <- flist[grepl(paste0("^", model, "$"), sub(".*_(iid|pi|rho|pirho)_.*", "\\1", flist))]
        # Read files
        setNames(lapply(model_files, function(file) {
            readRDS(file.path(data_path, file))
        }), paste0(sub(".*_([1-5]).Rds", "\\1", model_files)))
    }), model_vector)
}

alluvial_from_clustering <- function(clusterings) {
    bicl <- sapply(clusterings, function(clust) {
        sum(sapply(clust$partition, function(col) col$BICL))
    })
    bicl_df <- as.data.frame(bicl) |>
        rownames_to_column() |>
        mutate(axis = rowname)


    sapply(clusterings, function(clust) {
        clust$cluster
    }) |>
        as.data.frame() %>%
        rownames_to_column("id") %>% # Ajouter une colonne pour identifier les lignes
        pivot_longer(
            cols = as.character(seq(1, 5)), # Sélectionner dynamiquement les colonnes des clusters
            names_to = "axis", # Nom de la colonne pour les axes
            values_to = "cluster" # Nom de la colonne pour les valeurs
        ) -> clusters_long

    p_alluvial <- ggplot(data = clusters_long, aes(x = axis, stratum = cluster, alluvium = id)) +
        geom_flow(aes(fill = id)) +
        geom_stratum() +
        geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
        theme_minimal() +
        scale_fill_okabe_ito() +
        xlab("Clustering repetition")
    p_bicl <- ggplot(data = bicl_df, aes(x = axis, y = bicl)) +
        geom_point() +
        theme_minimal() +
        xlab("Clustering repetition")

    p <- p_alluvial + p_bicl +
        plot_layout(nrow = 2L) + plot_layout(guides = "collect")

    return(p)
}

alluvial_list_models <- function(clustering_list) {
    lapply(seq_along(clustering_list), function(i) {
        alluvial_from_clustering(clustering_list[[i]]) + ggtitle(paste0("Model: ", names(clustering_list)[i])) +
            plot_layout(tag_level = "new") + plot_annotation(tag_levels = c("1", "A"))
    })
}

baldock_list <- extract_from_flist(baldock_flist)
(baldock_clustering_plot <- wrap_plots(alluvial_list_models(baldock_list)))



traveset_flist <- flist[grepl("Traveset", flist)]
traveset_flist <- traveset_flist[order(factor(sub(".*_([1-5]).Rds", "\\1", traveset_flist), levels = seq(1, 5)))]
traveset_flist <- traveset_flist[order(factor(sub(".*_(iid|pi|rho|pirho)_.*", "\\1", traveset_flist), levels = model_vector))]
traveset_flist

traveset_list <- extract_from_flist(traveset_flist)
(traveset_clustering_plot <- wrap_plots(alluvial_list_models(traveset_list)) + plot_annotation(title = "Traveset"))


souza_flist <- flist[grepl("Souza", flist)]
souza_flist <- souza_flist[order(factor(sub(".*_([1-5]).Rds", "\\1", souza_flist), levels = seq(1, 5)))]
souza_flist <- souza_flist[order(factor(sub(".*_(iid|pi|rho|pirho)_.*", "\\1", souza_flist), levels = model_vector))]
souza_flist
souza_list <- extract_from_flist(souza_flist)
(souza_clustering_plot <- wrap_plots(alluvial_list_models(souza_list)) + plot_annotation(title = "Souza"))

cordeniz_flist <- flist[grepl("Cordeniz", flist)]
cordeniz_flist <- cordeniz_flist[order(factor(sub(".*_([1-5]).Rds", "\\1", cordeniz_flist), levels = seq(1, 5)))]
cordeniz_flist <- cordeniz_flist[order(factor(sub(".*_(iid|pi|rho|pirho)_.*", "\\1", cordeniz_flist), levels = model_vector))]
cordeniz_flist

cordeniz_list <- extract_from_flist(cordeniz_flist)
(cordeniz_clustering_plot <- wrap_plots(alluvial_list_models(cordeniz_list)) + plot_annotation(title = "Cordeniz"))

trojelsgaard_flist <- flist[grepl("Trojelsgaard", flist)]
trojelsgaard_flist <- trojelsgaard_flist[order(factor(sub(".*_([1-5]).Rds", "\\1", trojelsgaard_flist), levels = seq(1, 5)))]
trojelsgaard_flist <- trojelsgaard_flist[order(factor(sub(".*_(iid|pi|rho|pirho)_.*", "\\1", trojelsgaard_flist), levels = model_vector))]
trojelsgaard_flist
trojelsgaard_list <- extract_from_flist(trojelsgaard_flist)
(trojelsgaard_clustering_plot <- wrap_plots(alluvial_list_models(trojelsgaard_list)) + plot_annotation(title = "Trojelsgaard"))

gibson_flist <- flist[grepl("Gibson", flist)]
gibson_flist <- gibson_flist[order(factor(sub(".*_([1-5]).Rds", "\\1", gibson_flist), levels = seq(1, 5)))]
gibson_flist <- gibson_flist[order(factor(sub(".*_(iid|pi|rho|pirho)_.*", "\\1", gibson_flist), levels = model_vector))]
gibson_flist
gibson_list <- extract_from_flist(gibson_flist)
(gibson_clustering_plot <- wrap_plots(alluvial_list_models(gibson_list)) + plot_annotation(title = "Gibson"))


ggsave(here("figures", "applications", "subdore-per-author", "subdore_baldock_alluvial_clusterings.png"), baldock_clustering_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "subdore-per-author", "subdore_traveset_alluvial_clusterings.png"), traveset_clustering_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "subdore-per-author", "subdore_souza_alluvial_clusterings.png"), souza_clustering_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "subdore-per-author", "subdore_cordeniz_alluvial_clusterings.png"), cordeniz_clustering_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "subdore-per-author", "subdore_trojelsgaard_alluvial_clusterings.png"), trojelsgaard_clustering_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "subdore-per-author", "subdore_gibson_alluvial_clusterings.png"), gibson_clustering_plot, width = 12, height = 6, dpi = 300)


# Meso scale plots

wrap_plots(lapply(baldock_list[["iid"]][[1]]$partition, plot, type = "meso", values = TRUE)) + plot_annotation(title = "Baldock", subtitle = "iid")

baldock_list[["pi"]][[1]]$partition[[1]] |> plot(type = "meso", values = TRUE, mixture = TRUE)
