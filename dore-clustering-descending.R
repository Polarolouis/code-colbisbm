library(colSBM)
library(here)

set.seed(123L)
# Load dorebipartite data
data("dorebipartite")

# Saving all output to file
zz <- file(here("all.log"), open = "wt")
sink(zz)
sink(zz, type = "message")
# Clusterize
clustering <- clusterize_bipartite_networks(
    netlist = dorebipartite,
    colsbm_model = "iid",
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "no_mc"
    )
)
sink(type = "message")
sink()
file.show(here("all.log"))
