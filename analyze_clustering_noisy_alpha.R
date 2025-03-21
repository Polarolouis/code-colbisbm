library(here)
library(ggplot2)
library(ggokabeito)
library(tidyverse)

df <- readRDS(here("simulations", "clustering", "noisy_alpha", "noisy_alpha.Rds"))
noisy_alpha_plot <- df |>
    mutate(epsilon = as.factor(epsilon)) |>
    group_by(epsilon) |>
    select(epsilon, ari_truth) |>
    ggplot(aes(x = epsilon, y = ari_truth, group = epsilon, fill = epsilon)) +
    geom_boxplot() +
    labs(x = "$\\epsilon$", y = "ARI to truth") +
    # lims(y = c(0.5, 1)) +
    theme_minimal() +
    theme(aspect.ratio = 1L, axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0)) +
    ggtitle("ARI to truth by epsilon for noisy alpha")

ggsave(here("figures", "simulations", "clustering", "noisy_alpha", "noisy_alpha.png"), noisy_alpha_plot, width = 6, height = 6, dpi = 300)
