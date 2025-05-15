library(colSBM)
library(here)
library(stringr)
library(future.apply)
library(ggplot2)

baldock_matrices <- readRDS(here("data", "baldock2019-binary-matrices.Rds"))

fit <- estimate_colBiSBM(
    netlist = baldock_matrices,
    net_id = names(baldock_matrices),
    colsbm_model = "iid",
    global_opts = list(
        verbosity = 0,
        backend = "future"
    ),
    fit_opts = list(
        max_vem_steps = 10000L
    )
)
