library(colSBM)
library(gtools)
library(future.apply)
library(progressr)
set.seed(1234)
# Generate increasing group sizes
Q11 <- seq(1, 33, by = 1)
Q21 <- 3

rho1 <- c(0.2, 0.3, 0.5)

Q12 <- 3
Q22 <- 3

alpha2 <- matrix(c(
    0.9, 0.1, 0.1,
    0.2, 0.7, 0.1,
    0.1, 0.1, 0.8
), nrow = Q12)

pi2 <- c(0.2, 0.3, 0.5)
rho2 <- c(0.2, 0.3, 0.5)


#  Generate proportions matching Q11 and Q21 with a dirichlet distribution
pis1 <- lapply(Q11, function(Q1) {
    rdirichlet(1, alpha = rep(1, Q1))
})

alphas1 <- lapply(Q11, function(Q1) {
    matrix(runif(Q1 * Q21, 0.1, 1), nrow = Q1)
})

# time_taken <- sapply(seq_along(Q11), function(idx) {
#     start_time <- Sys.time()
#     graphon_distance_all_permutations(pis = list(pis1[[idx]], pi2), rhos = list(rho1, rho2), alphas = list(alphas1[[idx]], alpha2))
#     end_time <- Sys.time()
#     return(end_time - start_time)
# })
# Parallelized version
#  With a progressbar
with_progress({
    pb <- progressr::progressor(along = Q11)
    time_taken_parallel <- lapply(seq_along(Q11), function(idx) {
        start_time <- Sys.time()
        colSBM:::graphon_distance_marginals(pis = list(pis1[[idx]], pi2), rhos = list(rho1, rho2), alphas = list(matrix(alphas1[[idx]], nrow = idx), alpha2))
        pb()
        end_time <- Sys.time()
        return(end_time - start_time)
    })
})

library(ggplot2)

# Plot the time taken
time_taken_df <- data.frame(Q11 = Q11, time_taken = as.numeric(time_taken_parallel))
ggplot(time_taken_df, aes(x = Q11, y = time_taken)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(0, tail(Q11, 1), by = 1)) +
    labs(x = "Number of groups in network 1", y = "Time taken (s)") +
    theme_minimal()
