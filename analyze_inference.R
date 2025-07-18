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
    group_by(epsilon_alpha) |> # in grepl readd |un|ov
    select(which(!grepl("*_(ARI|BICL|secs)_*", colnames(averaged_print_data)),
        arr.ind = TRUE
    ))

if (!dir.exists(here("tables", "simulations", "inference"))) {
    dir.create(here("tables", "simulations", "inference"), recursive = TRUE)
}

length_col <- (ncol(averaged_print_data) - 1) / 4

all_names <- c(
    "$\\bm{1}_{\\widehat{Q_1} \\lt 4}$",
    "$\\bm{1}_{\\widehat{Q_1} = 4}$",
    "$\\bm{1}_{\\widehat{Q_1} \\gt 4}$",
    "$\\bm{1}_{\\widehat{Q_2} \\lt 4}$",
    "$\\bm{1}_{\\widehat{Q_2} = 4}$",
    "$\\bm{1}_{\\widehat{Q_2} \\gt 4}$"
)

eq_only_names <- c("$\\bm{1}_{\\widehat{Q_1} = 4}$", "$\\bm{1}_{\\widehat{Q_2} = 4}$")

(kbl(averaged_print_data,
    format = "html", booktabs = FALSE, escape = FALSE,
    col.names = c(
        "$\\epsilon_{\\alpha}$",
        rep(all_names, 4)
    ),
    linesep = "",
    vline = "|",
    # align = "|l|cc|cc|cc|cc|cccc|cccc|",
    caption = "The proportion of dataset where the correct number of blocks is selected.",
) |>
    kable_styling(font_size = 10L) |>
    add_header_above(c(" ", "iid" = length_col, "$\\\\pi$" = length_col, "$\\\\rho$" = length_col, "$\\\\pi\\\\rho$" = length_col), escape = FALSE, border_left = TRUE, border_right = TRUE, line_sep = 0)) |>
    save_kable(
        file = here(
            "tables", "simulations",
            "inference", "inference_table.tex"
        ),
        format = "latex",
        position = "H",
        size = "small",
        escape = FALSE
    )

## ----function_per_model, echo = FALSE------------------------------------------------------------------------------------------------------------------------------------------------
dataframe_per_model <- function(model) {
    averaged_print_data |>
        select(epsilon_alpha, starts_with(paste0(model, "_")))
}

## ----proportion-preferred_model, echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------


## ----proportion_preferred_figure, echo = FALSE---------------------------------------------------------------------------------------------------------------------------------------
#| fig.cap="\\label{fig:inference-proportion-preferred}Plot of the proportions of different preferred models in function of \\eps[\\alpha]",
#| fig.asp = 0.5,
#| fig.pos = "H",
#| fig.width = 7,
#| fig.height = 4,
#| dpi=300
output_tikz_folder <- here("tikz", "simulations", "inference")
if (!dir.exists(output_tikz_folder)) {
    dir.create(output_tikz_folder, recursive = TRUE)
}
proportion_preferred_data <- result_data_frame |>
    group_by(epsilon_alpha, preferred_model) |>
    summarise(n = n()) |>
    mutate(prop_model = n / sum(n)) |>
    ungroup() |>
    select(-n)

levels(proportion_preferred_data$preferred_model) <- c(
    "sep", "$iid$", "$\\pi$", "$\\rho$",
    "$\\pi\\rho$"
)
(plot <- proportion_preferred_data |>
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
    scale_color_okabe_ito(drop = FALSE) +
    scale_fill_okabe_ito(drop = FALSE) +
    xlab("$\\epsilon_{\\alpha}$") +
    ylab("Preferred model proportions") +
    theme_minimal() +
    theme(
        aspect.ratio = 1L,
        axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0),
        axis.text.y = element_text(size = 6)
    ) +
    geom_col(position = "stack"))

tikz(
    file = file.path(output_tikz_folder, "model-proportions.tex"), width = 4L,
    height = 3L,
    standAlone = TRUE
)

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
