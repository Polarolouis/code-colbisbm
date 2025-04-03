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


readRDS(here("data", "dore-binary-matrices.Rds")) -> all_dore_matrices
baldock_dore_matrices <- all_dore_matrices[grepl("Baldock", x = names(all_dore_matrices))]
plants_insects_lists <- lapply(baldock_dore_matrices, function(mat) {
    list(
        insects = rownames(mat),
        plants = colnames(mat)
    )
})

names(plants_insects_lists) <- names(plants_insects_lists) |> gsub(pattern = "Baldock201[1,9]_", replacement = "")

# Prepare latex table with merged common insects and plants
library(tibble)
library(kableExtra)
library(knitr)
library(purrr)
library(stringr)

outer(
    plants_insects_lists,
    plants_insects_lists,
    Vectorize(function(x, y) {
        paste0("$\\frac{", length(intersect(x$plants, y$plants)), "}{", length(union(x$plants, y$plants)), "},\\frac{", length(intersect(x$insects, y$insects)), "}{", length(union(x$insects, y$insects)), "}$")
    })
) -> nb_common

nb_common[upper.tri(nb_common, diag = TRUE)] <- ""
# Remove Baldock and keep cities from names
nb_common |>
    kable(format = "latex", caption = "Number of common insects/total insects and plants/total plants", escape = FALSE)
