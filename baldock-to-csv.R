library(stringr)
library(dplyr)

full_interaction_df <- read.table("data/dore/interaction-data.txt", header = TRUE, sep = "\t")

full_interaction_df |>
    filter(str_detect(web, pattern = "Baldock2019")) |>
    write.csv2("data/dore/baldock-interaction.csv", row.names = FALSE)
