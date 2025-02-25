library("colSBM")
library("sbm")
library("patchwork")
options(future.globals.maxSize = Inf)

set.seed(0L)

nr <- 120
nc <- 120

alpha <- matrix(c(
    0.8, 0.7, 0.25, 0.05, 0, 0,
    0.6, 0.4, 0.1, 0.05, 0, 0,
    0.4, 0.3, 0.9, 0.15, 0.05, 0.05,
    0.05, 0.05, 0.25, 0.8, 0.7, 0.25,
    0, 0, 0.05, 0.6, 0.4, 0.1,
    0, 0, 0.05, 0.4, 0.3, 0.9
), nrow = 6, byrow = TRUE)

alpha_up <- alpha[1:4, 1:4]
alpha_down <- alpha[3:6, 3:6]
alpha_center <- alpha[3:4, 3:4]

pi_ud <- c(0.1, 0.4, 0.2, 0.3)
pi_c <- c(0.3, 0.7)

rho_ud <- c(0.35, 0.1, 0.2, 0.35)
rho_c <- c(0.9, 0.1)

collection <- list(generate_bipartite_collection(
    nr, nc,
    pi_ud, rho_ud,
    alpha_up,
    M = 1,
    model = "iid",
    return_memberships = FALSE
)[[1]], generate_bipartite_collection(
    nr, nc,
    pi_c, rho_c,
    alpha_center,
    M = 1,
    model = "iid",
    return_memberships = FALSE
)[[1]], generate_bipartite_collection(
    nr, nc,
    pi_ud, rho_ud,
    alpha_down,
    M = 1,
    model = "iid",
    return_memberships = FALSE
)[[1]])

# Fitting sbm
fits_sbm <- lapply(collection, function(collection) {
    estimateBipartiteSBM(collection)
})
wrap_plots(lapply(fits_sbm, plot))

# Fitting colsbm
library("future.callr")
plan(tweak("callr", workers = 3L))
fit_colsbm_iid <- estimate_colBiSBM(
    netlist = collection,
    net_id = c("up", "center", "down"),
    colsbm_model = "iid"
)

fit_colsbm_pirho <- estimate_colBiSBM(
    netlist = collection,
    net_id = c("up", "center", "down"),
    colsbm_model = "pirho"
)

plot(fit_colsbm_iid$best_fit, type = "meso", mixture = TRUE)
