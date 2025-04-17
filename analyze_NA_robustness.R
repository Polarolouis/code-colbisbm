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
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")
# options(
#     tikzLatexPackages = c(getOption("tikzLatexPackages"), "\\usepackage{tikz}")
# )
filenames <- here("simulations", "NA_robustness", "NA_robustness_24-03-2025_13-04-01_1-800.Rds")

results_df <- readRDS(filenames)

matching_filenames <- filenames |>
    grep(pattern = "(modular|nested)", value = TRUE) |>
    head(n = 2)

results_df <- results_df |>
    select(-c(repetition))

auc_df <-
    results_df |>
    select(-contains(c("ari", "elapsed_secs"))) |>
    pivot_longer(
        cols = c(auc_LBM, auc_colBiSBM),
        values_to = "auc",
        names_prefix = "auc_",
        names_to = "method"
    ) |>
    mutate(model = ifelse(method == "LBM", paste0("sep-", model), model)) |>
    mutate_at(vars(model, struct), as.factor) |>
    mutate(model = forcats::fct_relevel(
        model,
        "iid", "sep-iid", "pi", "sep-pi", "rho", "sep-rho", "pirho", "sep-pirho"
    )) |>
    filter(prop_NAs <= 0.7)

levels(auc_df[["model"]]) <- levels(auc_df[["model"]]) |>
    str_replace(fixed("iid"), "$iid$") |>
    str_replace(fixed("pirho"), "$\\pi\\rho$") |>
    str_replace("pi$", "$\\\\pi$") |>
    str_replace("rho$", "$\\\\rho$")

library(patchwork)

lapply(c("iid", "\\pi", "\\rho", "\\pi\\rho"), function(model_name) {
    color <- switch(model_name,
        iid = 2,
        "\\pi" = 3,
        "\\rho" = 4,
        "\\pi\\rho" = 5
    )

    auc_df |>
        filter(str_equal(model, paste0("$", model_name, "$")) | str_equal(model, paste0("sep-$", model_name, "$"))) |>
        # Remove model_name from sep-model_name
        mutate(model = str_replace(model, fixed(paste0("sep-$", model_name, "$")), "sep")) |>
        ggplot() +
        aes(x = factor(prop_NAs), y = auc, fill = model, color = model) +
        geom_boxplot(notch = TRUE) +
        labs(x = "$p_{\\texttt{NA}}$", y = "AUC") +
        scale_fill_okabe_ito(order = c(color, 1), alpha = 0.5) +
        scale_color_okabe_ito(order = c(color, 1)) +
        guides(
            fill = guide_legend(title = "Model"),
            color = guide_legend(title = "Model")
        ) +
        ylim(c(0.5, 1)) +
        theme_minimal()
}) -> auc_plots

(auc_plot <- wrap_plots(auc_plots) + plot_layout(axis_titles = "collect"))

output_tikz_folder <- here(
    "tikz", "simulations",
    "na_robustness"
)
if (!dir.exists(output_tikz_folder)) {
    dir.create(output_tikz_folder, recursive = TRUE)
}

tikz(
    file = file.path(output_tikz_folder, "auc-model.tex"), width = 6.5,
    height = 4.5,
    standAlone = TRUE
)
print(auc_plot)
dev.off()
## ARI plots
ari_df <- results_df |>
    select(-contains(c("auc", "elapsed_secs"))) |>
    group_by(prop_NAs, model) |>
    pivot_longer(
        cols = c(arirow_LBM:aricol_colBiSBM),
        values_to = "ari", names_to = c("dim", "method"),
        names_pattern = "(row|col)_(LBM|colBiSBM)"
    ) |>
    mutate(model = ifelse(method == "LBM", paste0("sep-", model), model)) |>
    mutate_at(vars(model, struct, dim, method), as.factor) |>
    mutate(model = forcats::fct_relevel(
        model,
        "iid", "sep-iid", "pi", "sep-pi", "rho", "sep-rho", "pirho", "sep-pirho"
    )) |>
    filter(prop_NAs <= 0.7)

levels(ari_df[["model"]]) <- levels(ari_df[["model"]]) |>
    str_replace(fixed("iid"), "$iid$") |>
    str_replace(fixed("pirho"), "$\\pi\\rho$") |>
    str_replace("pi$", "$\\\\pi$") |>
    str_replace("rho$", "$\\\\rho$")

ari_df[["dim"]] <- relevel(ari_df[["dim"]], ref = "row")
dim.labs <- c("Rows : $d = 1$", "Columns : $d = 2$")
names(dim.labs) <- c("row", "col")

lapply(c("iid", "\\pi", "\\rho", "\\pi\\rho"), function(model_name) {
    color <- switch(model_name,
        iid = 2,
        "\\pi" = 3,
        "\\rho" = 4,
        "\\pi\\rho" = 5
    )

    ari_df |>
        filter(str_equal(model, paste0("$", model_name, "$")) | str_equal(model, paste0("sep-$", model_name, "$"))) |>
        # Remove model_name from sep-model_name
        mutate(model = str_replace(model, fixed(paste0("sep-$", model_name, "$")), "sep")) |>
        ggplot() +
        aes(x = factor(prop_NAs), y = ari, fill = model, color = model) +
        geom_boxplot(notch = TRUE) +
        labs(x = "$p_{\\texttt{NA}}$", y = "ARI") +
        guides(
            fill = guide_legend(title = "Model"),
            color = guide_legend(title = "Model")
        ) +
        scale_fill_okabe_ito(order = c(color, 1), alpha = 0.5) +
        scale_color_okabe_ito(order = c(color, 1)) +
        # facet_grid(struct ~ dim, labeller = labeller(dim = dim.labs)) +
        theme_minimal() +
        theme(axis.text.x = element_text(size = 5))
}) -> ari_plots
(ari_plot <- wrap_plots(ari_plots) + plot_layout(axis_titles = "collect"))
tikz(
    file = file.path(output_tikz_folder, "ari-dim-model.tex"), width = 6.5,
    height = 4,
    standAlone = TRUE
)
print(ari_plot)
dev.off()
