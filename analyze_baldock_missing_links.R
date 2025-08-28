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
filenames <- here(data_dir, "missing_links_baldock_1752848977.14268.Rds")

results_df <- readRDS(filenames)

results_df <- rbind(results_df, readRDS(here(data_dir, "missing_links_sep.Rds")))

results_df <- results_df |>
    select(-c(repetitions))

auc_df <-
    results_df |>
    mutate_at(vars(model), as.factor) |>
    mutate(model = forcats::fct_relevel(
        model,
        "sep", "iid", "pi", "rho", "pirho"
    )) |>
    mutate(missing_replacement = ifelse(is.na(missing_replacement), "Missing Dyads", "Missing Links"))

levels(auc_df[["model"]]) <- levels(auc_df[["model"]]) |>
    str_replace(fixed("iid"), "$iid$") |>
    str_replace(fixed("pirho"), "$\\pi\\rho$") |>
    str_replace("pi$", "$\\\\pi$") |>
    str_replace("rho$", "$\\\\rho$")

(auc_plot <- ggplot(auc_df) +
    aes(y = AUC, x = as.factor(epsilon), fill = model) +
    geom_boxplot(aes(fill = model)) +
    facet_wrap(. ~ missing_replacement, ncol = 2L) +
    scale_fill_okabe_ito() +
    labs(fill = "Model") +
    xlab("$p_{\\mbox{mis}}$") +
    ylab("ROC AUC") +
    ylim(0.8, 0.91) +
    theme_bw() +
    scale_y_continuous(n.breaks = 5) +
    theme(aspect.ratio = 0.8))

fig_folder <- here(
    "figures", "simulations",
    "missing_links"
)
if (!dir.exists(fig_folder)) {
    dir.create(fig_folder, recursive = TRUE)
}

tikz(
    file = file.path(fig_folder, "auc-model.tex"), width = 6.5,
    height = 4.5,
    standAlone = TRUE
)
print(auc_plot)
dev.off()

png(filename = file.path(fig_folder, "auc-model.png"), width = 1200)
print(auc_plot)
dev.off()
