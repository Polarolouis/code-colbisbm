library("colSBM")
library("ggplot2")
library("tikzDevice")
library("patchwork")

my_meso_plot <- function(best_fit) {
    plot(best_fit, type = "meso", values = TRUE, mixture = TRUE, text_size = 5.5) + plot_layout(guides = "collect") & theme(legend.position = "bottom", text = element_text(family = "serif"), legend.title = element_text(size = 12), axis.text = element_text(size = 11))
}


baldock_clustering <- readRDS("data/baldock_clustering.Rds")
baldock_clustering$partition[[1]]$net_id <- "Kenyan-TB+JN"

(kenyan_iid_plot <- my_meso_plot(baldock_clustering$partition[[1]]) & guides(fill = "none"))

options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")
tikz(file = "figures/applications/baldock/baldock-kenyan-iid.tex", width = 10.3, height = 7.53, standAlone = TRUE)
kenyan_iid_plot
dev.off()
