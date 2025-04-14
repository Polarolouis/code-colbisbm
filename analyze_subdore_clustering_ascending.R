library(here)
library(ggplot2)
library(patchwork)

models <- c("iid")

datalist <- readRDS(here("applications", "subdore-ascending", "subdore-ascending.Rds"))
lapply(datalist, function(clustering) {
    clustering$partition <- lapply(clustering$partition, function(col) {
        class(col) <- c("fitBipartiteSBMPop", "R6")
        col
    })
    clustering
}) -> datalist

nb_rep <- length(datalist) / length(models)

names(datalist) <- paste0(rep(models), ".", rep(seq(1, nb_rep), each = length(models)))

iid_clust_log <- grepl("iid", names(datalist))
clusterings_iid <- datalist[iid_clust_log]

library(dplyr)
library(tidyr)
library(tibble)
library("ggalluvial")
library("ggokabeito")



sapply(clusterings_iid, function(clust) {
    clust$cluster
}) |>
    as.data.frame() %>%
    mutate(authors = ifelse(rownames(.) |> grepl(pattern = "Baldock", fixed = TRUE), "Baldock", ifelse(rownames(.) |> grepl(pattern = "Gibson", fixed = TRUE), "Gibson", ifelse(rownames(.) |> grepl(pattern = "Souza", fixed = TRUE), "Souza", ifelse(rownames(.) |> grepl(pattern = "Trojelsgaard", fixed = TRUE), "Trojelsgaard", ifelse(rownames(.) |> grepl(pattern = "Cordeniz", fixed = TRUE), "Cordeniz", ifelse(rownames(.) |> grepl(pattern = "Traveset", fixed = TRUE), "Traveset", "Other"))))))) -> iid_clusters

iid_bicl <- sapply(clusterings_iid, function(clust) {
    tail(clust$bicl_history, 1)
})
iid_bicl_df <- as.data.frame(iid_bicl) |> rownames_to_column()


as.data.frame(iid_clusters) %>%
    rownames_to_column("id") %>% # Ajouter une colonne pour identifier les lignes
    pivot_longer(
        cols = starts_with("iid."), # Sélectionner dynamiquement les colonnes des clusters
        names_to = "axis", # Nom de la colonne pour les axes
        values_to = "cluster" # Nom de la colonne pour les valeurs
    ) -> iid_clusters_long

iid_clusters_long <- left_join(iid_clusters_long, iid_bicl_df, by = c("axis" = "rowname")) %>%
    mutate(bicl = iid_bicl) %>%
    select(-iid_bicl)


ggplot(data = iid_clusters_long, aes(x = axis, stratum = cluster, alluvium = id, fill = authors)) +
    geom_flow() +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_okabe_ito() +
    theme_minimal() +
    labs(
        title = "Alluvial diagram for sub-Doré networks",
        subtitle = "iid"
    ) -> iid_alluvial

ggplot(iid_clusters_long, aes(x = axis, y = bicl)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> iid_bicl_plot
iid_alluvial / iid_bicl_plot -> iid_alluvial
iid_alluvial
