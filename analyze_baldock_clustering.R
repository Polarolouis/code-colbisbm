library(here)
library(ggplot2)
library(patchwork)

datalist <- readRDS(here("applications", "baldock", "baldock.Rds"))

iid_clust_log <- grepl("iid", names(datalist))
clusterings_iid <- datalist[iid_clust_log]

plot_all_meso <- function(clusterings) {
    wrap_plots(lapply(clusterings, function(clust) {
        wrap_plots(
            lapply(clust$partition, function(col) {
                plot(col, type = "meso", values = TRUE) + theme(legend.position = "none")
            }),
        )
    })) + plot_layout(tag_level = "new") + plot_annotation(tag_levels = c("1", "A"))
}
plot_all_meso(clusterings_iid) + plot_layout(nrow = 3L) + plot_annotation(title = "Meso-scale partition for Baldock networks", subtitle = "iid") -> iid_plot

pi_clust_log <- grepl("pi", names(datalist))
clusterings_pi <- datalist[pi_clust_log]

plot_all_meso(clusterings_pi) + plot_annotation(title = "Meso-scale partition for Baldock networks", subtitle = "$\\pi$") -> pi_plot

rho_clust_log <- grepl("rho", names(datalist))
clusterings_rho <- datalist[rho_clust_log]

plot_all_meso(clusterings_rho) + plot_annotation(title = "Meso-scale partition for Baldock networks", subtitle = "$\\rho$") -> rho_plot

pirho_clust_log <- grepl("pirho", names(datalist))
clusterings_pirho <- datalist[pirho_clust_log]

pirho_plot <- plot_all_meso(clusterings_pirho) + plot_annotation(title = "Meso-scale partition for Baldock networks", subtitle = "$\\pi\\rho$")

ggsave(here("figures", "applications", "baldock", "baldock_meso_iid.png"), iid_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "baldock", "baldock_meso_pi.png"), pi_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "baldock", "baldock_meso_rho.png"), rho_plot, width = 12, height = 6, dpi = 300)
ggsave(here("figures", "applications", "baldock", "baldock_meso_pirho.png"), pirho_plot, width = 12, height = 6, dpi = 300)
