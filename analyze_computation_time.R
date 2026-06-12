save_path <- "simulations/computation_time"

rds_files <- list.files(save_path, include.dirs = FALSE, pattern = "Rds")

computation_complex_time_df <- readRDS(file.path(save_path, rds_files[1]))
library(ggplot2)

library(tidyverse)
computation_complex_time_df %>%
    group_by(Q1, Q2) %>%
    select(-rep) %>%
    summarise_at(.vars = c("duration"), list(avg = mean, sd = sd)) %>%
    ggplot(aes(x = Q1, y = Q2, fill = avg)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(
        x = "Q1",
        y = "Q2",
        fill = "Time"
    ) +
    theme_minimal()
