library(colSBM)
library(future.apply)
library(future.callr)
plan(callr, workers = 3L)

set.seed(123L)

nr <- 75
nc <- 75


alpha <- matrix(c(
    0.8, 0.5, 0.2, 0.05,
    0.45, 0.2, 0.05, 0.05,
    0.35, 0.05, 0.8, 0.05,
    0.05, 0.05, 0.05, 0.65
), nrow = 4L, byrow = TRUE)

alpha1 <- alpha[1:3, 1:3]
alpha2 <- alpha[2:4, 2:4]

pi <- c(0.45, 0.3, 0.25)
rho <- c(0.45, 0.2, 0.35)

netlist <- c(generate_bipartite_collection(
    nr = nr, nc = nc,
    pi = pi, rho = rho,
    M = 1, model = "iid",
    alpha = alpha1
), generate_bipartite_collection(
    nr = nr, nc = nc,
    pi = pi, rho = rho,
    M = 1, model = "iid",
    alpha = alpha2
))

names(netlist) <- c("alpha1", "alpha2")

fit <- estimate_colBiSBM(
    netlist = netlist,
    colsbm_model = "pirho",
    net_id = names(netlist),
    distribution = "bernoulli",
    fit_opts = list(
        algo_ve = "fp",
        minibatch = TRUE,
        verbosity = 0
    ),
    global_opts = list(
        verbosity = 2,
        backend = "future"
    ),
    nb_run = 3L
)
