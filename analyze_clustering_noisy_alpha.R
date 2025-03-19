library(here)
library(ggplot2)
library(ggokabeito)
library(tidyverse)

df <- readRDS(here("simulations", "clustering", "noisy_alpha", "noisy_alpha.Rds"))
df |>
    mutate(epsilon = as.factor(epsilon)) |>
    group_by(epsilon) |>
    select(epsilon, ari_truth) |>
    ggplot(aes(x = epsilon, y = ari_truth, group = epsilon, fill = epsilon)) +
    geom_boxplot() +
    labs(x = expression(epsilon), y = "ARI to truth") +
    theme_minimal() +
    theme(aspect.ratio = 1L, axis.text.x = element_text(angle = -45, vjust = .5, hjust = 0)) +
    ggtitle("ARI to truth by epsilon for noisy alpha")
