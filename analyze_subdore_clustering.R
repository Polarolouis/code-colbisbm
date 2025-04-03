library(here)
library(ggplot2)
library(patchwork)

datalist <- readRDS(here("applications", "subdore", "subdore.Rds"))
lapply(datalist, function(clustering) {
    clustering$partition <- lapply(clustering$partition, function(col) {
        class(col) <- c("fitBipartiteSBMPop", "R6")
        col
    })
    clustering
}) -> datalist

nb_rep <- length(datalist) / 4L

names(datalist) <- paste0(rep(c("iid", "pi", "rho", "pirho")), ".", rep(seq(1, nb_rep), each = 4))

iid_clust_log <- grepl("iid", names(datalist))
clusterings_iid <- datalist[iid_clust_log]

plot_all_meso <- function(clusterings) {
    wrap_plots(lapply(clusterings, function(clust) {
        wrap_plots(
            lapply(clust$partition, function(col) {
                plot(col, type = "meso", values = TRUE) + theme(legend.position = "none")
            }),
        ) + plot_layout(tag_level = "new")
    })) + plot_annotation(tag_levels = c("A", "1"))
}
iid_plot <- plot_all_meso(clusterings_iid) + plot_layout(nrow = 3L) + plot_annotation(title = "Meso-scale partition for sub-Doré networks", subtitle = "iid")

pi_clust_log <- grepl("pi.[0-9]{,2}$", names(datalist))
clusterings_pi <- datalist[pi_clust_log]

plot_all_meso(clusterings_pi) + plot_annotation(title = "Meso-scale partition for sub-Doré networks", subtitle = "$\\pi$") -> pi_plot

rho_clust_log <- grepl("rho.[0-9]{,2}$", names(datalist))
clusterings_rho <- datalist[rho_clust_log]

plot_all_meso(clusterings_rho) + plot_annotation(title = "Meso-scale partition for sub-Doré networks", subtitle = "$\\rho$") -> rho_plot

pirho_clust_log <- grepl("pirho.[0-9]{,2}$", names(datalist))
clusterings_pirho <- datalist[pirho_clust_log]

pirho_plot <- plot_all_meso(clusterings_pirho) + plot_annotation(title = "Meso-scale partition for sub-Doré networks", subtitle = "$\\pi\\rho$")

# ggsave(here("figures", "applications", "subdore", "subdore_meso_iid.png"), iid_plot, width = 12, height = 6, dpi = 300)
# ggsave(here("figures", "applications", "subdore", "subdore_meso_pi.png"), pi_plot, width = 12, height = 6, dpi = 300)
# ggsave(here("figures", "applications", "subdore", "subdore_meso_rho.png"), rho_plot, width = 12, height = 6, dpi = 300)
# ggsave(here("figures", "applications", "subdore", "subdore_meso_pirho.png"), pirho_plot, width = 12, height = 6, dpi = 300)

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
    sum(sapply(clust$partition, function(col) col$BICL))
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

sapply(clusterings_pi, function(clust) {
    clust$cluster
}) |>
    as.data.frame() %>%
    mutate(authors = ifelse(rownames(.) |> grepl(pattern = "Baldock", fixed = TRUE), "Baldock", ifelse(rownames(.) |> grepl(pattern = "Gibson", fixed = TRUE), "Gibson", ifelse(rownames(.) |> grepl(pattern = "Souza", fixed = TRUE), "Souza", ifelse(rownames(.) |> grepl(pattern = "Trojelsgaard", fixed = TRUE), "Trojelsgaard", ifelse(rownames(.) |> grepl(pattern = "Cordeniz", fixed = TRUE), "Cordeniz", ifelse(rownames(.) |> grepl(pattern = "Traveset", fixed = TRUE), "Traveset", "Other"))))))) -> pi_clusters

pi_bicl <- sapply(clusterings_pi, function(clust) {
    sum(sapply(clust$partition, function(col) col$BICL))
})
pi_bicl_df <- as.data.frame(pi_bicl) |> rownames_to_column()

as.data.frame(pi_clusters) %>%
    rownames_to_column("id") %>% # Ajouter une colonne pour identifier les lignes
    pivot_longer(
        cols = starts_with("pi."), # Sélectionner dynamiquement les colonnes des clusters
        names_to = "axis", # Nom de la colonne pour les axes
        values_to = "cluster" # Nom de la colonne pour les valeurs
    ) -> pi_clusters_long

pi_clusters_long <- left_join(pi_clusters_long, pi_bicl_df, by = c("axis" = "rowname")) %>%
    mutate(bicl = pi_bicl) %>%
    select(-pi_bicl)

ggplot(data = pi_clusters_long, aes(x = axis, stratum = cluster, alluvium = id, fill = authors)) +
    geom_flow() +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_okabe_ito() +
    theme_minimal() +
    labs(
        title = "Alluvial diagram for sub-Doré networks",
        subtitle = "$\\pi$"
    ) -> pi_alluvial

ggplot(pi_clusters_long, aes(x = axis, y = bicl)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> pi_bicl_plot
pi_alluvial / pi_bicl_plot -> pi_alluvial
pi_alluvial

sapply(clusterings_rho, function(clust) {
    clust$cluster
}) |>
    as.data.frame() %>%
    mutate(authors = ifelse(rownames(.) |> grepl(pattern = "Baldock", fixed = TRUE), "Baldock", ifelse(rownames(.) |> grepl(pattern = "Gibson", fixed = TRUE), "Gibson", ifelse(rownames(.) |> grepl(pattern = "Souza", fixed = TRUE), "Souza", ifelse(rownames(.) |> grepl(pattern = "Trojelsgaard", fixed = TRUE), "Trojelsgaard", ifelse(rownames(.) |> grepl(pattern = "Cordeniz", fixed = TRUE), "Cordeniz", ifelse(rownames(.) |> grepl(pattern = "Traveset", fixed = TRUE), "Traveset", "Other"))))))) -> rho_clusters

rho_bicl <- sapply(clusterings_rho, function(clust) {
    sum(sapply(clust$partition, function(col) col$BICL))
})
rho_bicl_df <- as.data.frame(rho_bicl) |> rownames_to_column()

as.data.frame(rho_clusters) %>%
    rownames_to_column("id") %>% # Ajouter une colonne pour identifier les lignes
    pivot_longer(
        cols = starts_with("rho."), # Sélectionner dynamiquement les colonnes des clusters
        names_to = "axis", # Nom de la colonne pour les axes
        values_to = "cluster" # Nom de la colonne pour les valeurs
    ) -> rho_clusters_long

rho_clusters_long <- left_join(rho_clusters_long, rho_bicl_df, by = c("axis" = "rowname")) %>%
    mutate(bicl = rho_bicl) %>%
    select(-rho_bicl)

ggplot(data = rho_clusters_long, aes(x = axis, stratum = cluster, alluvium = id, fill = authors)) +
    geom_flow() +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_okabe_ito() +
    theme_minimal() +
    labs(
        title = "Alluvial diagram for sub-Doré networks",
        subtitle = "$\\rho$"
    ) -> rho_alluvial

ggplot(rho_clusters_long, aes(x = axis, y = bicl)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> rho_bicl_plot
rho_alluvial / rho_bicl_plot -> rho_alluvial



sapply(clusterings_pirho, function(clust) {
    clust$cluster
}) |>
    as.data.frame() %>%
    mutate(authors = ifelse(rownames(.) |> grepl(pattern = "Baldock", fixed = TRUE), "Baldock", ifelse(rownames(.) |> grepl(pattern = "Gibson", fixed = TRUE), "Gibson", ifelse(rownames(.) |> grepl(pattern = "Souza", fixed = TRUE), "Souza", ifelse(rownames(.) |> grepl(pattern = "Trojelsgaard", fixed = TRUE), "Trojelsgaard", ifelse(rownames(.) |> grepl(pattern = "Cordeniz", fixed = TRUE), "Cordeniz", ifelse(rownames(.) |> grepl(pattern = "Traveset", fixed = TRUE), "Traveset", "Other"))))))) -> pirho_clusters

pirho_bicl <- sapply(clusterings_pirho, function(clust) {
    sum(sapply(clust$partition, function(col) col$BICL))
})
pirho_bicl_df <- as.data.frame(pirho_bicl) |> rownames_to_column()

as.data.frame(pirho_clusters) %>%
    rownames_to_column("id") %>% # Ajouter une colonne pour identifier les lignes
    pivot_longer(
        cols = starts_with("pirho."), # Sélectionner dynamiquement les colonnes des clusters
        names_to = "axis", # Nom de la colonne pour les axes
        values_to = "cluster" # Nom de la colonne pour les valeurs
    ) -> pirho_clusters_long

pirho_clusters_long <- left_join(pirho_clusters_long, pirho_bicl_df, by = c("axis" = "rowname")) %>%
    mutate(bicl = pirho_bicl) %>%
    select(-pirho_bicl)


ggplot(data = pirho_clusters_long, aes(x = axis, stratum = cluster, alluvium = id, fill = authors)) +
    geom_flow() +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_okabe_ito() +
    theme_minimal() +
    # Rotate the x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        title = "Alluvial diagram for sub-Doré networks",
        subtitle = "$\\pi\\rho$"
    ) -> pirho_alluvial
ggplot(pirho_clusters_long, aes(x = axis, y = bicl)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> pirho_bicl_plot

pirho_alluvial / pirho_bicl_plot -> pirho_alluvial



iid_alluvial
pi_alluvial
rho_alluvial
pirho_alluvial
