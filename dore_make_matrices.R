library(here)
library(dplyr)

dore_df <- read.table(here("data", "dore", "interaction-data.txt"), sep = "\t", header = TRUE)

dore_df <- dore_df |> mutate(plantfull = paste(plantorder, plantfamily, plantgenus, plantspecies, sep = "_"), insectfull = paste(insectorder, insectfamily, insectgenus, insectspecies, sep = "_"))

# Matrices at network level
dore_networks_matrices <- setNames(lapply(unique(dore_df$id_network), function(id) {
    adj_table <- table(dore_df[dore_df$id_network == id, c("insectfull", "plantfull")])
    return(matrix(adj_table, ncol = ncol(adj_table), dimnames = dimnames(adj_table)))
}), nm = unique(dore_df$web))

saveRDS(dore_networks_matrices, here("data", "dore", "dore_networks_matrices.Rds"))

# Matrices at aggregated network level
aggreg_names <- sapply(unique(dore_df$id_network_aggreg), function(id) {
    names <- unique(dore_df[dore_df$id_network_aggreg == id, "web"])
    merged_names <- ifelse(length(names) > 1L, paste(names, collapse = "-"), names)
})
dore_aggregated_networks_matrices <- setNames(lapply(unique(dore_df$id_network_aggreg), function(id) {
    adj_table <- table(dore_df[dore_df$id_network_aggreg == id, c("insectfull", "plantfull")])
    return(matrix(adj_table, ncol = ncol(adj_table), dimnames = dimnames(adj_table)))
}), nm = aggreg_names)

saveRDS(dore_aggregated_networks_matrices, here("data", "dore", "dore_aggregated_networks_matrices.Rds"))
