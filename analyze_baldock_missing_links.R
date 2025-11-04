## ----libraries, echo = FALSE, include = FALSE----------------------------------------------------------------------------------------------------------------------------------------
library("ggplot2")
library("ggokabeito")
library("tidyr")
library("dplyr")
library("stringr")
library("knitr")
library("kableExtra")
library("vctrs")
library("stringr")
library("here")
library("tikzDevice")
library("forcats")
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")
# options(
#     tikzLatexPackages = c(getOption("tikzLatexPackages"), "\\usepackage{tikz}")
# )
data_dir <- here("simulations", "missing_links")
# filenames <- here(data_dir, "missing_links_baldock_1752848977.14268.Rds")
filenames <- here(data_dir, "missing_links_baldock_1761742524.6736.Rds")


results_df <- readRDS(filenames)

results_df <- rbind(results_df, readRDS(here(data_dir, "missing_links_sep.Rds")))

results_df <- results_df |>
    select(-c(repetitions))

vgae_df <- read.csv("data/dore/baldock_vgae_metrics_per_epsilon.csv")



vgae_df <- vgae_df |>
    rename(AUC = AUC, epsilon = Epsilon, possible_missing_network = Missing.Network) |>
    mutate(model = "VGAE", missing_replacement = 0) |>
    select(names(results_df))

ggplot(vgae_df) +
    aes(x = as.factor(epsilon), y = AUC, color = as.factor(possible_missing_network)) +
    geom_boxplot()

sep_vgae_df <- read.csv("data/dore/baldock_sep_vgae_metrics_per_epsilon.csv")
sep_vgae_df <- sep_vgae_df |>
    rename(AUC = AUC, epsilon = Epsilon, possible_missing_network = Missing.Network) |>
    mutate(model = "sepVGAE", missing_replacement = 0) |>
    select(names(results_df))
ggplot(sep_vgae_df) +
    aes(x = as.factor(epsilon), y = AUC, color = as.factor(possible_missing_network)) +
    geom_boxplot()


results_df <- rbind(results_df, vgae_df, sep_vgae_df)

auc_df <-
    results_df |>
    mutate_at(vars(model), as.factor) |>
    mutate(model = forcats::fct_relevel(
        model,
        "sep", "iid", "pi", "rho", "pirho"
    )) |>
    mutate(missing_replacement = ifelse(is.na(missing_replacement), "Missing Dyads", "Missing Links"))

levels(auc_df[["model"]]) <- levels(auc_df[["model"]]) |>
    str_replace(fixed("iid"), "iid") |>
    str_replace(fixed("pirho"), "$\\pi\\rho$") |>
    str_replace("pi$", "$\\\\pi$") |>
    str_replace("rho$", "$\\\\rho$")

(auc_plot <- ggplot(auc_df) +
    aes(y = AUC, x = as.factor(epsilon), fill = model, color = model) +
    geom_boxplot(aes(fill = model), notch = TRUE) +
    facet_wrap(. ~ missing_replacement,
        ncol = 2L,
        scales = "free"
    ) +
    scale_color_okabe_ito() +
    scale_fill_okabe_ito(alpha = 0.5) +
    labs(fill = "Model", color = "Model") +
    xlab("$p_{\\mbox{mis}}$") +
    ylab("ROC AUC") +
    # ylim(0.8, 0.91) +
    theme_minimal() +
    scale_y_continuous(n.breaks = 5) +
    theme(legend.position = "bottom"))

(auc_plot_missing_links <- ggplot(auc_df |> filter(missing_replacement == "Missing Links")) +
    aes(y = AUC, x = as.factor(epsilon), fill = model, color = model) +
    geom_boxplot(aes(fill = model), notch = TRUE) +
    scale_color_okabe_ito() +
    scale_fill_okabe_ito(alpha = 0.5) +
    labs(fill = "Model", color = "Model") +
    xlab("$p_{\\mbox{mis}}$") +
    ylab("ROC AUC") +
    ylim(0.8, 0.91) +
    theme_minimal() +
    scale_y_continuous(n.breaks = 5) +
    theme(legend.position = "bottom"))

fig_folder <- here(
    "figures", "simulations",
    "missing_links"
)
if (!dir.exists(fig_folder)) {
    dir.create(fig_folder, recursive = TRUE)
}

tikz(
    file = file.path(fig_folder, "auc-model.tex"), width = 6.5,
    height = 3,
    standAlone = TRUE
)
print(auc_plot)
dev.off()

tikz(
    file = file.path(fig_folder, "auc-model-missing-links.tex"), width = 10,
    height = 6.5,
    standAlone = TRUE
)
print(auc_plot_missing_links)
dev.off()

png(filename = file.path(fig_folder, "auc-model.png"), width = 1200)
print(auc_plot)
dev.off()
