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

pis_clust_log <- grepl("pi.[0-9]{,1}$", names(datalist))
clusterings_pi <- datalist[pi_clust_log]

plot_all_meso(clusterings_pi) + plot_annotation(title = "Meso-scale partition for sub-Doré networks", subtitle = "$\\pi$") -> pi_plot

rho_clust_log <- grepl("rho.[0-9]{,1}$", names(datalist))
clusterings_rho <- datalist[rho_clust_log]

plot_all_meso(clusterings_rho) + plot_annotation(title = "Meso-scale partition for sub-Doré networks", subtitle = "$\\rho$") -> rho_plot

pirho_clust_log <- grepl("pirho.[0-9]{,1}$", names(datalist))
clusterings_pirho <- datalist[pirho_clust_log]

pirho_plot <- plot_all_meso(clusterings_pirho) + plot_annotation(title = "Meso-scale partition for sub-Doré networks", subtitle = "$\\pi\\rho$")

# ggsave(here("figures", "applications", "subdore", "subdore_meso_iid.png"), iid_plot, width = 12, height = 6, dpi = 300)
# ggsave(here("figures", "applications", "subdore", "subdore_meso_pi.png"), pi_plot, width = 12, height = 6, dpi = 300)
# ggsave(here("figures", "applications", "subdore", "subdore_meso_rho.png"), rho_plot, width = 12, height = 6, dpi = 300)
# ggsave(here("figures", "applications", "subdore", "subdore_meso_pirho.png"), pirho_plot, width = 12, height = 6, dpi = 300)

sapply(clusterings_iid, function(clust) {
    clust$cluster
}) -> iid_clusters

sapply(clusterings_pi, function(clust) {
    clust$cluster
}) -> pi_clusters

sapply(clusterings_rho, function(clust) {
    clust$cluster
}) -> rho_clusters

sapply(clusterings_pirho, function(clust) {
    clust$cluster
}) -> pirho_clusters
library("ggalluvial")



ggplot(data = as.data.frame(iid_clusters)) +
    aes(axis1 = iid.1, axis2 = iid.2, axis3 = iid.3, axis4 = iid.4, axis5 = iid.5) +
    geom_alluvium(aes(fill = as.factor(iid.1))) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Alluvial diagram for sub-Doré networks", subtitle = "iid") -> iid_alluvial
ggplot(data = as.data.frame(pi_clusters)) +
    aes(axis1 = pi.1, axis2 = pi.2, axis3 = pi.3, axis4 = pi.4, axis5 = pi.5) +
    geom_alluvium(aes(fill = as.factor(pi.1))) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Alluvial diagram for sub-Doré networks", subtitle = "$\\pi$") -> pi_alluvial
ggplot(data = as.data.frame(rho_clusters)) +
    aes(axis1 = rho.1, axis2 = rho.2, axis3 = rho.3, axis4 = rho.4, axis5 = rho.5) +
    geom_alluvium(aes(fill = as.factor(rho.1))) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Alluvial diagram for sub-Doré networks", subtitle = "$\\rho$") -> rho_alluvial

ggplot(data = as.data.frame(pirho_clusters)) +
    aes(axis1 = pirho.1, axis2 = pirho.2, axis3 = pirho.3, axis4 = pirho.4, axis5 = pirho.5) +
    geom_alluvium(aes(fill = as.factor(pirho.1))) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Alluvial diagram for sub-Doré networks", subtitle = "$\\pi\\rho$") -> pirho_alluvial

iid_alluvial + pi_alluvial + rho_alluvial + pirho_alluvial -> alluvial_plot
