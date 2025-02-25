## ----libraries, echo = FALSE, include = FALSE----------------------------------------------------------------------------------------------------------------------------------------
library("ggplot2")
library("ggokabeito")
library("tidyr")
library("dplyr")
library("stringr")
library("knitr")
library("kableExtra")
library("stringr")
library("here")
library("latex2exp")
library("tikzDevice")
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")



## ----impoting-data, echo = FALSE----------------------------------------------------------------------------------------------------------------------------
# data_list <- lapply(filenames, function(file) lapply(readRDS(file), function(model) model$list_clustering))
df_netclust <- readRDS("simulations/clustering/9collection/descending.Rds")
df_netclust$clust <- "descending"
df_ascending <- readRDS("simulations/clustering/9collection/ascending.Rds")
df_ascending$clust <- "ascending"
df_iid_on_pirho <- readRDS("simulations/clustering/9collection/iid_on_pirho.Rds")
df_iid_on_pirho$clust <- "descending"

df_netclust <- rbind(df_netclust, df_iid_on_pirho)

df_netclust$model <- factor(df_netclust$model, levels = c(
    "iid", "pi",
    "rho", "pirho",
    "iid_on_pirho"
))

df_netclust$clust <- factor(df_netclust$clust, levels = c(
    "ascending", "descending"
))

df_netclust$model <- df_netclust$model |>
    case_match("iid" ~ "$iid$",
        "pi" ~ "$\\pi$",
        "rho" ~ "$\\rho$",
        "pirho" ~ "$\\pi\\rho$",
        "iid_on_pirho" ~ "$iid$ on $\\pi\\rho$",
        .ptype = factor(levels = c(
            "$iid$", "$\\pi$",
            "$\\rho$", "$\\pi\\rho$",
            "$iid$ on $\\pi\\rho$"
        ))
    )

## ----netclustering-ARI-boxplot, echo = FALSE----------------------------------------------------------------------------------------------------------------
#| dpi = 300,
#| fig.asp = 0.5,
#| fig.cap = "\\label{}ARI of the partition obtained by clustering in function of $\\eps$"
ari_plot <- ggplot(df_netclust) +
    aes(x = as.factor(epsilon), y = ARI) +
    scale_color_okabe_ito(order = c(1L, 8L)) +
    scale_fill_okabe_ito(order = 2L:9L, alpha = 0.5) +
    xlab("$\\epsilon_{\\alpha}$") +
    guides(
        fill = guide_legend(title = "Model"),
        color = guide_legend(title = "Clustering")
    ) +
    ylab("ARI of the clustering") +
    geom_boxplot(aes(fill = model, color = clust), notch = TRUE) +
    theme_minimal()

bicl_plot <- ggplot(df_netclust) +
    aes(x = as.factor(epsilon), y = BICL) +
    scale_color_okabe_ito(order = c(1L, 8L)) +
    scale_fill_okabe_ito(order = 2L:9L, alpha = 0.5) +
    xlab("$\\epsilon_{\\alpha}$") +
    guides(
        fill = guide_legend(title = "Model"),
        color = guide_legend(title = "Clustering")
    ) +
    ylab("BICL of the clustering") +
    geom_boxplot(aes(fill = model, color = clust)) +
    theme_minimal()

output_tikz_folder <- here(
    "figures",
    "simulations", "clustering", "9collection"
)

if (!dir.exists(output_tikz_folder)) {
    dir.create(output_tikz_folder, recursive = TRUE)
}

# tikz(
#     file = file.path(output_tikz_folder, "ari-clustering.tex"), width = 5L,
#     height = 2,
#     standAlone = TRUE
# )
# print(ari_plot)
# dev.off()

png(
    file = file.path(output_tikz_folder, "ari-clustering-descending.png"),
    width = 960,
    height = 960
)
print(ari_plot)
dev.off()

png(
    file = file.path(output_tikz_folder, "bicl-clustering-descending.png"),
    width = 960,
    height = 960
)
print(bicl_plot)
dev.off()
