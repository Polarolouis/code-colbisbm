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
library("tikzDevice")
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")


meanse <- function(x, ...) {
    mean1 <- signif(round(mean(x, na.rm = T), 2), 5) # calculate mean and round
    se1 <- signif(round(sd(x, na.rm = T) / sqrt(sum(!is.na(x))), 2), 2) # std error - round adding zeros
    out <- paste(mean1, "$\\pm$", se1) # paste together mean plus/minus and standard error
    if (str_detect(out, "NA")) {
        out <- "NA"
    } # if missing do not add plusminus
    if (se1 == 0) {
        out <- paste(mean1)
    }
    return(out)
}


## ----import-data, echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------
# filenames <- list.files(
#     path = here("code", "results", "simulations", "inference", "bernoulli"),
#     # pattern = "",
#     full.names = TRUE
# )

filename <- here("simulations", "inference", "bernoulli_inference.Rds")

col_id_BICLS <- c(11, 16, 23, 30, 37)
result_data_frame <- readRDS(filename)

# Compute the preferred model
result_data_frame <- cbind(result_data_frame, preferred_model = sapply(seq_len(nrow(result_data_frame)), function(n) sub("_BICL", "", names(which.max(result_data_frame[n, col_id_BICLS])))))

result_data_frame$preferred_model <- factor(result_data_frame$preferred_model, levels = c(
    "sep", "iid", "pi",
    "rho", "pirho"
))

result_data_frame <- result_data_frame |>
    mutate(
        iid_Q1_un = as.integer(iid_Q1 < 4L),
        iid_Q1_eq = as.integer(iid_Q1 == 4L),
        iid_Q1_ov = as.integer(iid_Q1 > 4L),
        iid_Q2_un = as.integer(iid_Q2 < 4L),
        iid_Q2_eq = as.integer(iid_Q2 == 4L),
        iid_Q2_ov = as.integer(iid_Q2 > 4L),
        pi_Q1_un = as.integer(pi_Q1 < 4L),
        pi_Q1_eq = as.integer(pi_Q1 == 4L),
        pi_Q1_ov = as.integer(pi_Q1 > 4L),
        pi_Q2_un = as.integer(pi_Q2 < 4L),
        pi_Q2_eq = as.integer(pi_Q2 == 4L),
        pi_Q2_ov = as.integer(pi_Q2 > 4L),
        rho_Q1_un = as.integer(rho_Q1 < 4L),
        rho_Q1_eq = as.integer(rho_Q1 == 4L),
        rho_Q1_ov = as.integer(rho_Q1 > 4L),
        rho_Q2_un = as.integer(rho_Q2 < 4L),
        rho_Q2_eq = as.integer(rho_Q2 == 4L),
        rho_Q2_ov = as.integer(rho_Q2 > 4L),
        pirho_Q1_un = as.integer(pirho_Q1 < 4L),
        pirho_Q1_eq = as.integer(pirho_Q1 == 4L),
        pirho_Q1_ov = as.integer(pirho_Q1 > 4L),
        pirho_Q2_un = as.integer(pirho_Q2 < 4L),
        pirho_Q2_eq = as.integer(pirho_Q2 == 4L),
        pirho_Q2_ov = as.integer(pirho_Q2 > 4L)
    ) |>
    select(-c(
        iid_Q1, iid_Q2, pi_Q1, pi_Q2,
        rho_Q1, rho_Q2, pirho_Q1, pirho_Q2
    ))

# TODO Finir d'ajouter les colonnes sur le blocs


## ----inference_table, echo = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
averaged_print_data <- result_data_frame |>
    group_by(epsilon_alpha) |>
    summarise(across(-preferred_model, list("avrg" = meanse))) |>
    select(-c(2:10))
averaged_print_data <- averaged_print_data |>
    select(which(!grepl("*_BICL_*", colnames(averaged_print_data)),
        arr.ind = TRUE
    ))

## ----function_per_model, echo = FALSE------------------------------------------------------------------------------------------------------------------------------------------------
dataframe_per_model <- function(model) {
    averaged_print_data |>
        select(epsilon_alpha, starts_with(paste0(model, "_")))
}

## ----per_model_table, echo = FALSE, results='asis', message=FALSE, warning = FALSE---------------------------------------------------------------------------------------------------

dir.create(here(
    "mia-rapport-2024", "tables", "simulations",
    "inference"
), recursive = TRUE)

for (model in c("sep", "iid", "pi", "rho", "pirho")) {
    kable_ari_colnames <- c(
        "$\\eps[\\alpha]$", # "BIC-L",
        "$\\overline{\\text{ARI}}_{1}$",
        "$\\overline{\\text{ARI}}_{2}$", "$\\text{ARI}_{1}$", "$\\text{ARI}_{2}$"
    )
    model_name <- model
    if (model != "sep") {
        kable_blocrecov_colnames <- c(
            "$\\eps[\\alpha]$",
            "$\\mathbbb{1}_{\\widehat{Q_1}<Q_1}$",
            "$\\mathbbb{1}_{\\widehat{Q_1}=Q_1}$",
            "$\\mathbbb{1}_{\\widehat{Q_1}>Q_1}$",
            "$\\mathbbb{1}_{\\widehat{Q_2}<Q_2}$",
            "$\\mathbbb{1}_{\\widehat{Q_2}=Q_2}$",
            "$\\mathbbb{1}_{\\widehat{Q_2}>Q_2}$"
        )
    }
    if (model == "pirho") {
        model_name <- "$\\pi\\rho$"
    } else {
        if (model != "iid" && model != "sep") {
            model_name <- paste0("$\\", model, "$")
        } else {
            model_name <- paste0("$", model, "$")
        }
    }
    full_dataframe <- dataframe_per_model(model)
    ari_dataframe <- full_dataframe |>
        group_by(epsilon_alpha) |>
        select(-contains(c("Q1", "Q2")))
    blocrecov_dataframe <- full_dataframe |>
        group_by(epsilon_alpha) |>
        select(contains(c("Q1", "Q2")))

    ari_kable <- kable(ari_dataframe,
        escape = FALSE,
        booktabs = TRUE,
        digits = 2,
        caption = paste0(
            "\\label{subtab:ari_per_model_", model,
            "}Quality metrics for ",
            ifelse(model != "sep", paste0(model_name, "$\\text{-colBiSBM}$"), "$sep\\text{-BiSBM}$")
        ),
        col.names = kable_ari_colnames,
        format = "latex"
    ) |>
        kable_styling(latex_options = "scale_down")
    if (model != "sep") {
        blocrecov_kable <- kable(blocrecov_dataframe,
            escape = FALSE,
            booktabs = TRUE,
            digits = 2,
            caption = paste0(
                "\\label{subtab:blocrecov_per_model_", model,
                "}Bloc recovery for ",
                ifelse(model != "sep", paste0(model_name, "$\\text{-colBiSBM}$"), "$sep\\text{-BiSBM}$")
            ),
            col.names = kable_blocrecov_colnames,
            format = "latex"
        ) |>
            kable_styling(latex_options = "scale_down")
        both_kables <- kables(list(ari_kable, blocrecov_kable))
    } else {
        both_kables <- ari_kable
    }
    both_kables <- both_kables |>
        gsub("\\begin{table}", "\\begin{subtable}{\\textwidth}", x = _, fixed = TRUE) |>
        gsub("\\end{table}", "\\end{subtable}", x = _, fixed = TRUE)
    cat("",
        "\\begin{table}[H]",
        "\\centering",
        paste0("\\caption{\\label{tab:inference_results_", model, "}Inference results for ", model_name, "}"),
        both_kables,
        "\\end{table}",
        "",
        sep = "\n",
        file = here(
            "mia-rapport-2024", "tables", "simulations",
            "inference", paste0(model, ".tex")
        )
    )
}


## ----proportion-preferred_model, echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------
proportion_preferred_data <- result_data_frame |>
    group_by(epsilon_alpha, preferred_model) |>
    summarise(n = n()) |>
    mutate(prop_model = n / sum(n)) |>
    ungroup() |>
    select(-n)

proportion_preferred_table <- proportion_preferred_data |>
    pivot_wider(
        names_from = preferred_model,
        values_from = prop_model, values_fill = 0
    )

cat(kable(proportion_preferred_table,
    escape = FALSE,
    booktabs = TRUE,
    digits = 2,
    position = "!h",
    caption = "\\label{tab:proportion-preferred-table}Proportions of models selected per \\eps[\\alpha] (data for Figure \\ref{fig:inference-proportion-preferred})",
    col.names = c(
        "\\eps[\\alpha]",
        "$sep\\text{-}BiSBM$",
        "$iid\\text{-}colBiSBM$",
        "$\\pi\\text{-}colBiSBM$",
        "$\\rho\\text{-}colBiSBM$",
        "$\\pi\\rho\\text{-}colBiSBM$"
    ),
    align = "rccccc",
    format = "latex"
) |>
    kable_styling(latex_options = "scale_down"), file = here(
    "mia-rapport-2024", "tables", "simulations",
    "inference", "preferred.tex"
))


## ----proportion_preferred_figure, echo = FALSE---------------------------------------------------------------------------------------------------------------------------------------
#| fig.cap="\\label{fig:inference-proportion-preferred}Plot of the proportions of different preferred models in function of \\eps[\\alpha]",
#| fig.asp = 0.5,
#| fig.pos = "H",
#| fig.width = 7,
#| fig.height = 4,
#| dpi=300
output_tikz_folder <- here("mia-rapport-2024", "tikz", "simulations", "inference")
if (!dir.exists(output_tikz_folder)) {
    dir.create(output_tikz_folder, recursive = TRUE)
}

tikz(
    file = file.path(output_tikz_folder, "model-proportions.tex"), width = 4L,
    height = 3L,
    standAlone = TRUE
)
levels(proportion_preferred_data$preferred_model) <- c(
    "sep", "$iid$", "$\\pi$", "$\\rho$",
    "$\\pi\\rho$"
)
plot <- proportion_preferred_data |>
    ggplot() +
    aes(
        x = epsilon_alpha, y = prop_model, color = preferred_model,
        fill = preferred_model
    ) +
    guides(
        fill = guide_legend(title = "Model"),
        color = guide_legend(title = "Model")
    ) +
    scale_x_continuous(breaks = seq(from = 0.0, to = 0.24, by = 0.03)) +
    scale_color_okabe_ito() +
    scale_fill_okabe_ito() +
    xlab("$\\epsilon_{\\alpha}$") +
    ylab("Preferred model proportions") +
    theme_minimal() +
    theme(
        aspect.ratio = 1L,
        axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0),
        axis.text.y = element_text(size = 6)
    ) +
    geom_col(position = "stack")
print(plot)
dev.off()

averaged_data <- result_data_frame |>
    group_by(epsilon_alpha) |>
    select(-contains(c("BICL", "pi1", "rho2", "Q1", "Q2", "elapsed_secs", "preferred_model", "repetition"))) |>
    pivot_longer(contains("ARI"), names_pattern = "([a-zA-Z\\_]*)_ARI", values_to = "ARI") |>
    separate(name, sep = "_", into = c("model", "ARI_type", "dim")) |>
    mutate(
        model = forcats::fct_relevel(
            model,
            "sep", "iid", "pi", "rho", "pirho"
        ),
        ARI_type = forcats::fct_relevel(ARI_type, "mean", "double"),
        dim = forcats::fct_relevel(dim, "row", "col")
    ) |>
    ungroup() |>
    group_by(epsilon_alpha, model, ARI_type, dim) |>
    summarise(across(everything(), list("mean" = mean, "sd" = sd)))

unaveraged_data <- result_data_frame |>
    group_by(epsilon_alpha) |>
    select(-contains(c("BICL", "pi1", "rho2", "Q1", "Q2", "elapsed_secs", "preferred_model", "repetition"))) |>
    pivot_longer(contains("ARI"), names_pattern = "([a-zA-Z\\_]*)_ARI", values_to = "ARI") |>
    separate(name, sep = "_", into = c("model", "ARI_type", "dim")) |>
    mutate(
        model = forcats::fct_relevel(
            model,
            "sep", "iid", "pi", "rho", "pirho"
        ),
        ARI_type = forcats::fct_relevel(ARI_type, "mean", "double"),
        dim = forcats::fct_relevel(dim, "row", "col")
    ) |>
    ungroup() |>
    group_by(epsilon_alpha, model, ARI_type, dim)

dim.labs <- c("$d = 1$", "$d = 2$")
names(dim.labs) <- c("row", "col")
ARI_type.labs <- c("$\\overline{\\mbox{ARI}}_d$", "$\\mbox{ARI}_d$")
names(ARI_type.labs) <- c("mean", "double")

levels(averaged_data$model) <- c(
    "sep", "$iid$", "$\\pi$", "$\\rho$",
    "$\\pi\\rho$"
)
levels(unaveraged_data$model) <- c(
    "sep", "$iid$", "$\\pi$", "$\\rho$",
    "$\\pi\\rho$"
)

ARI_plots <- ggplot(averaged_data) +
    aes(x = epsilon_alpha, y = ARI_mean, color = model) +
    geom_point(aes(fill = model)) +
    geom_line() +
    geom_ribbon(
        aes(
            ymin = ARI_mean - ARI_sd,
            ymax = ARI_mean + ARI_sd, fill = model
        ),
        alpha = 0.05
    ) +
    scale_color_okabe_ito() +
    scale_fill_okabe_ito() +
    facet_grid(ARI_type ~ dim, labeller = labeller(
        ARI_type = ARI_type.labs,
        dim = dim.labs
    )) +
    guides(
        fill = guide_legend(title = "Model"),
        color = guide_legend(title = "Model")
    ) +
    labs(y = "", x = "$\\epsilon_{\\alpha}$") +
    theme_minimal() +
    theme(aspect.ratio = 1L, axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0))

# ggplot(unaveraged_data) +
#     aes(x = as.factor(epsilon_alpha), y = ARI, color = model) +
#     geom_boxplot()+
#     scale_color_okabe_ito() +
#     scale_fill_okabe_ito() +
#     facet_grid(ARI_type ~ dim, labeller = labeller(
#         ARI_type = ARI_type.labs,
#         dim = dim.labs
#     )) +
#     guides(
#         fill = guide_legend(title = "Model"),
#         color = guide_legend(title = "Model")
#     ) +
#     labs(y = "", x = "$\\epsilon_{\\alpha}$") +
#     theme_minimal() +
#     theme(aspect.ratio = 1L, axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0))

tikz(
    file = file.path(output_tikz_folder, "ari-plots.tex"), width = 6L,
    height = 5L,
    standAlone = TRUE
)
print(ARI_plots)
dev.off()
