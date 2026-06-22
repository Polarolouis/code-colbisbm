library(here)
library(tidyverse)
library(latex2exp)
save_path <- "simulations/computation_time"

M <- 10
nr <- nc <- 300

pattern <- paste0("M_", M, "_nr_", nr, "_nc_", nc)

rds_files <- list.files(save_path, include.dirs = FALSE, pattern = pattern)

computation_complex_time_df <- do.call("rbind", lapply(file.path(save_path, rds_files), readRDS)) %>% mutate(model = factor(model, levels = c("iid", "pirho"), labels = c("iid", "pirho")))

# computation_complex_time_df <- readRDS(file.path(save_path, rds_files[1]))
library(ggplot2)
library(ggrepel)
library(tidyverse)


avg_computation <- computation_complex_time_df %>%
    group_by(Q1, Q2, model) %>%
    select(-rep) %>%
    summarise_at(.vars = c("duration"), list(avg = function(x) as.numeric(mean(x)), sd = sd)) %>%
    mutate(labelQ2 = if_else(Q1 == max(Q1), as.character(Q2), NA_character_), labelQ1 = if_else(Q2 == max(Q2), as.character(Q1), NA_character_))
ggplot(avg_computation, aes(x = Q1, y = Q2, fill = avg)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "red") +
    labs(
        x = "Q1",
        y = "Q2",
        fill = "Time"
    ) +
    facet_wrap(~model) +
    theme_minimal()

model_colors <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2")
library(latex2exp)
(q1_comp_plot <- ggplot(computation_complex_time_df %>% filter(Q2 == 8), aes(x = Q1, fill = model)) +
    geom_boxplot(aes(y = duration, group = interaction(Q1, model)), notch = TRUE) +
    geom_line(data = avg_computation %>% filter(Q2 == 8), aes(y = avg, color = model)) +
    scale_fill_manual(values = model_colors[c(1, 4)], drop = FALSE, labels = c("iid", TeX("$\\pi\\rho$-colBiSBM"))) +
    scale_color_manual(values = model_colors[c(1, 4)], drop = FALSE, labels = c("iid", TeX("$\\pi\\rho$-colBiSBM"))) +
    scale_x_continuous(breaks = seq(2, 8)) +
    scale_y_continuous(n.breaks = 10) +
    labs(fill = "Model", color = "Model", y = "Time (s)", x = TeX("$Q_1$")) +
    ggtitle(TeX(paste0("$n_1=n_2=", nr, ",M=", M, "$"))) +
    theme_minimal())

(q2_comp_plot <- ggplot(computation_complex_time_df %>% filter(Q1 == 8), aes(x = Q2, fill = model)) +
    geom_boxplot(aes(y = duration, group = interaction(Q2, model)), notch = TRUE) +
    geom_line(data = avg_computation %>% filter(Q1 == 8), aes(y = avg, color = model)) +
    scale_fill_manual(values = model_colors[c(1, 4)], drop = FALSE, labels = c("iid", TeX("$\\pi\\rho$-colBiSBM"))) +
    scale_color_manual(values = model_colors[c(1, 4)], drop = FALSE, labels = c("iid", TeX("$\\pi\\rho$-colBiSBM"))) +
    scale_x_continuous(breaks = seq(2, 8)) +
    scale_y_continuous(n.breaks = 10) +
    labs(fill = "Model", color = "Model", y = "Time (s)", x = TeX("$Q_2$")) +
    theme_minimal())

library(patchwork)
library(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")
output_tikz_folder <- here(
    "tikz", "simulations",
    "computation_time"
)
if (!dir.exists(output_tikz_folder)) {
    dir.create(output_tikz_folder, recursive = TRUE)
}

width_tikz <- 10
height_tikz <- 3

pdf(
    file = file.path(output_tikz_folder, paste0("computation-time-", pattern, ".pdf")), width = width_tikz, height = height_tikz
)
q1_comp_plot + q2_comp_plot + plot_layout(guides = "collect", axis_titles = "collect")
dev.off()
