library(here)
library(ggplot2)
library(ggokabeito)
library(tidyverse)

df <- readRDS(here("simulations", "clustering", "noisy_links", "noisy_links.Rds"))
noisy_links_plot <- df |>
    mutate(epsilon = as.factor(epsilon)) |>
    group_by(epsilon) |>
    filter(type == "noisy") |>
    select(epsilon, ari_truth) |>
    ggplot(aes(x = epsilon, y = ari_truth, group = epsilon, fill = epsilon)) +
    geom_boxplot() +
    labs(x = "$\\epsilon$", y = "ARI to truth") +
    theme_minimal() +
    theme(aspect.ratio = 1L, axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0)) +
    ggtitle("ARI to truth by epsilon for noisy links")

clear_links_plot <- df |>
    mutate(epsilon = as.factor(epsilon)) |>
    group_by(epsilon) |>
    filter(type == "clear") |>
    select(epsilon, ari_truth) |>
    ggplot(aes(x = epsilon, y = ari_truth, group = epsilon, fill = epsilon)) +
    geom_boxplot() +
    labs(x = "$\\epsilon$", y = "ARI to truth") +
    theme_minimal() +
    theme(aspect.ratio = 1L, axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0)) +
    ggtitle("ARI to truth by epsilon for clear data")


ggsave(here("figures", "simulations", "clustering", "noisy_links", "noisy_links.png"), noisy_links_plot, width = 6, height = 6, dpi = 300)

ggsave(here("figures", "simulations", "clustering", "noisy_links", "clear_links.png"), clear_links_plot, width = 6, height = 6, dpi = 300)
