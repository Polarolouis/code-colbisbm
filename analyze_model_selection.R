## ----libraries, echo = FALSE, include = FALSE-----------------------------------------------------------------------------------------------------------------------
require("ggplot2")
require("ggokabeito")
require("knitr")
require("kableExtra")
require("stringr")
require("tidyr")
require("dplyr")
require("here")
require("ggh4x")
require("tikzDevice")
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")


## ----setup, echo = FALSE, include= FALSE----------------------------------------------------------------------------------------------------------------------------
options(knitr.table.format = "latex")

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


## ----import-data, echo = FALSE--------------------------------------------------------------------------------------------------------------------------------------
filename <- here("simulations", "model_selection", "model_selection.Rds")
result_data_frame <- readRDS(filename)
result_data_frame$preferred_model <- factor(result_data_frame$preferred_model, levels = c(
    "iid", "pi",
    "rho", "pirho"
))

# Adding a column accounting for true model iid, pi, rho or pirho
# result_data_frame <- result_data_frame %>% mutate(true_model = if (all( c(epsilon_pi >0, epsilon_rho > 0) == c(TRUE, TRUE))) print("pirho") else if (all( c(epsilon_pi >0, epsilon_rho > 0) == c(TRUE, FALSE))) print("pi") else if (all( c(epsilon_pi >0, epsilon_rho > 0) == c(F, T))) print("rho") else print("iid"))



## ----compute-table, echo = FALSE, include = FALSE-------------------------------------------------------------------------------------------------------------------
result_data_frame <- result_data_frame %>%
    mutate(
        iid_Q1 = as.integer(iid_Q1 == 3L),
        iid_Q2 = as.integer(iid_Q2 == 3L),
        pi_Q1 = as.integer(pi_Q1 == 3L),
        pi_Q2 = as.integer(pi_Q2 == 3L),
        rho_Q1 = as.integer(rho_Q1 == 3L),
        rho_Q2 = as.integer(rho_Q2 == 3L),
        pirho_Q1 = as.integer(pirho_Q1 == 3L),
        pirho_Q2 = as.integer(pirho_Q2 == 3L)
    )

model_comparison_eps_pi_rho <- result_data_frame %>%
    group_by(epsilon_pi, epsilon_rho, preferred_model) %>%
    summarise(n = n()) %>%
    mutate(prop_model = n / sum(n)) %>%
    select(-n)

model_comparison_eps_pi <- result_data_frame %>%
    group_by(epsilon_pi, preferred_model) %>%
    summarise(n = n(), rec_Q1 = mean(iid_Q1 + pi_Q1 + rho_Q1 + pirho_Q1) / 4) %>%
    mutate(prop_model = n / sum(n))

model_comparison_eps_rho <- result_data_frame %>%
    group_by(epsilon_rho, preferred_model) %>%
    summarise(n = n(), rec_Q2 = mean(iid_Q2 + pi_Q2 + rho_Q2 + pirho_Q2) / 4) %>%
    mutate(prop_model = n / sum(n))

bloc_recovery_df <- result_data_frame |>
    filter(epsilon_rho %in% c(0, 0.14, 0.28)) |>
    group_by(epsilon_pi, epsilon_rho) |>
    select(-(contains(c("BICL", "pi2", "rho2", "pi3", "rho3", "repetition", "elapsed_secs", "preferred_model")))) |>
    summarise_all(meanse)
# relocate(c(iid_Q2, pi_Q2, rho_Q2), .after = pirho_Q1)


model_proportion_df <- result_data_frame |>
    group_by(epsilon_pi, epsilon_rho, preferred_model) |>
    select(-(contains(c("BICL", "pi", "rho", "iid", "repetition", "elapsed_secs")))) |>
    summarise(n = n()) |>
    mutate(prop_model = n / sum(n)) |>
    select(-n) |>
    pivot_wider(
        names_from = preferred_model, values_from = prop_model, names_prefix = "proportion_" # ,
        # values_fill = 0
    )

table_df <- full_join(x = bloc_recovery_df, y = model_proportion_df)

## ----epsilon_plot, echo = FALSE, include = FALSE--------------------------------------------------------------------------------------------------------------------
#| fig.asp = 0.5,
#| fig.pos = "H",
#| fig.width = 7,
#| fig.height = 4,
#| dpi=300
levels(model_comparison_eps_pi$preferred_model) <- c(
    "sep", "$iid$", "$\\pi$", "$\\rho$",
    "$\\pi\\rho$"
)
levels(model_comparison_eps_rho$preferred_model) <- c(
    "sep", "$iid$", "$\\pi$", "$\\rho$",
    "$\\pi\\rho$"
)
levels(model_comparison_eps_pi_rho$preferred_model) <- c(
    "$iid$", "$\\pi$", "$\\rho$",
    "$\\pi\\rho$"
)


plot_pi_rho <- ggplot(model_comparison_eps_pi_rho, aes(
    x = "", y = prop_model,
    fill = preferred_model
)) +
    geom_col(position = "stack") +
    labs(
        title = "",
        y = "Selection proportions",
        x = "",
        fill = "Model"
    ) +
    scale_y_continuous(breaks = c(0, 1)) +
    scale_color_okabe_ito(order = 2L:9L) +
    scale_fill_okabe_ito(order = 2L:9L) +
    facet_nested("$\\epsilon_{\\rho}$" + factor(epsilon_rho, levels = rev(unique(epsilon_rho))) ~ "$\\epsilon_{\\pi}$" + epsilon_pi) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6), legend.position = "bottom")
plot_pi_rho
output_tikz_folder <- here(
    "tikz", "simulations",
    "model_selection"
)
if (!dir.exists(output_tikz_folder)) {
    dir.create(output_tikz_folder, recursive = TRUE)
}

tikz(
    file = file.path(output_tikz_folder, "eps-pi-rho-preferred.tex"), width = 6,
    height = 4,
    standAlone = TRUE
)
print(plot_pi_rho)
dev.off()

# table_folder <- here(
#     "mia-rapport-2024", "tables", "simulations",
#     "model_selection"
# )

# if (!dir.exists(table_folder)) {
#     dir.create(table_folder)
# }

# format_values <- function(x) {
#     ifelse(x == 1, "1", sprintf("%.3f", x))
# }

# print_table_df <- as.data.frame(lapply(table_df, format_values)) |>
#     rowwise() |>
#     mutate(across(proportion_iid:proportion_pirho, ~ ifelse(. == max(c_across(proportion_iid:proportion_pirho), na.rm = TRUE), paste0("\\textbf{", ., "}"), .))) %>%
#     ungroup() |>
#     rowwise() |>
#     mutate(across(iid_Q1:pirho_Q2, ~ ifelse(. == min(c_across(iid_Q1:pirho_Q2), na.rm = TRUE) & length(unique(c_across(iid_Q1:pirho_Q2))) > 1L, paste0("\\textbf{", ., "}"), .))) %>%
#     ungroup()



# options(knitr.kable.NA = "")
# cat(
#     print_table_df |>
#         # filter((epsilon_rho < 0.11 | epsilon_rho > 0.21)) |>
#         # filter((epsilon_pi < 0.11 | epsilon_pi > 0.21)) |>
#         # filter(!(epsilon_pi %in% c(0.11, 0.14, 0.18, 0.21))) |>
#         kable(
#             format = "latex", col.names = c(
#                 "$\\epsilon_{\\pi}$", "$\\epsilon_{\\rho}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_1}_{iid}=3}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_2}_{iid}=3}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_1}_{\\pi}=3}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_2}_{\\pi}=3}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_1}_{\\rho}=3}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_2}_{\\rho}=3}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_1}_{\\pi\\rho}=3}$",
#                 "$\\mathbbb{1}_{\\widehat{Q_2}_{\\pi\\rho}=3}$",
#                 "$iid$",
#                 "$\\pi$",
#                 "$\\rho$",
#                 "$\\pi\\rho$"
#             ), escape = FALSE,
#             caption = "\\label{tab:model-selection}Filtered block recovery and model selection proportions",
#             booktabs = TRUE, longtable = TRUE,
#             digits = 3L, full_width = TRUE
#         ) |>
#         column_spec(11L, border_left = TRUE) |>
#         add_header_above(c(
#             " " = 2L, "$iid$" = 2L,
#             "$\\\\pi$" = 2L, "$\\\\rho$" = 2L,
#             "$\\\\pi\\\\rho$" = 2L, " " = 4L
#         ), escape = FALSE) |>
#         add_header_above(c(
#             " " = 2L, "Block number recovery" = 8L,
#             "Model selection proportions" = 4L
#         )) |>
#         kable_styling(latex_options = c("repeat_header")) |>
#         collapse_rows(columns = 1:2, valign = "middle"),
#     file = file.path(table_folder, "model-selection.tex")
# )

if (!dir.exists(here("tables", "simulations", "model_selection"))) {
    dir.create(here("tables", "simulations", "model_selection"), recursive = TRUE)
}

kbl(bloc_recovery_df,
    format = "latex",
    escape = FALSE,
    linesep = "",
    borders = "",
    col.names = c(
        "$\\epsilon_{\\pi}$", "$\\epsilon_{\\rho}$",
        "$\\bm{1}_{\\widehat{Q_1}_{iid}=3}$",
        "$\\bm{1}_{\\widehat{Q_2}_{iid}=3}$",
        "$\\bm{1}_{\\widehat{Q_1}_{\\pi}=3}$",
        "$\\bm{1}_{\\widehat{Q_2}_{\\pi}=3}$",
        "$\\bm{1}_{\\widehat{Q_1}_{\\rho}=3}$",
        "$\\bm{1}_{\\widehat{Q_2}_{\\rho}=3}$",
        "$\\bm{1}_{\\widehat{Q_1}_{\\pi\\rho}=3}$",
        "$\\bm{1}_{\\widehat{Q_2}_{\\pi\\rho}=3}$"
    ),
    align = "|ll|cc|cc|cc|cc|c|",
    caption = "\\textit{From fixed  to varying block proportions}. Proportion of dataset where the correct number of blocks is
        selected. Only shows results for $\\epsilon_{\\rho}\\in \\{0, 0.14, 0.28\\}$.\\label{tab:model-selection-block-recovery}",
    position = "!ht",
) |>
    add_header_above(c(
        " " = 2L, "$iid$" = 2L,
        "$\\\\pi$" = 2L, "$\\\\rho$" = 2L,
        "$\\\\pi\\\\rho$" = 2L
    ), escape = FALSE, border_left = TRUE, border_right = TRUE) |>
    kable_styling(font_size = 9) |>
    save_kable(
        file = here("tables", "simulations", "model_selection", "block_recovery.tex"),
        format = "latex",
        latex_options = c("repeat_header")
    )
