library("colSBM")
library("future")
library("future.callr")

plan("callr")

nb_rep <- 3
max_Q <- 8
M <- 3

nr <- 300
nc <- 300

seq_Q <- seq(max_Q)

Q1 <- 5
Q <- c(Q1, Q1)
alpha_iid <- diag(exp(-1 / seq(Q1)))
alpha_pirho <- diag(exp(-1 / seq(Q1))) + abs(rnorm(n = Q1^2, sd = 0.01))

pis <- exp(-1 / seq(Q1))
pis <- pis / sum(pis)
rhos <- exp(-1 / rev(seq(Q1)))
rhos <- rhos / sum(rhos)

check_identifiability_iid(nr = rep(nr, M), nc = rep(nc, M), Q = Q, alpha = alpha_iid, pi = pis, rho = rhos)

collection_params <- generate_bipartite_collection(nr = nr, nc = nc, pi = pis, rho = rhos, alpha = alpha_pirho, M = M, model = "pirho", return_blockProp = TRUE)

collection <- collection_params[["collection"]]
pis <- collection_params[["blockProps"]][["pi"]]
rhos <- collection_params[["blockProps"]][["rho"]]


check_identifiability_pirho(nr = sapply(collection, nrow), nc = sapply(collection, ncol), Q = Q, alpha = alpha_pirho, pis = pis, rhos = rhos)

datasets <- lapply(seq(nb_rep), function(idx) {
    generate_bipartite_collection()
})
