library(colSBM)
library(gtools)

# Generate multiple parameters sets of fixed size
# To test if the order of graphon distance from all_permutations is the same as the order of graphon_distance_marginals
set.seed(1234)
Q11 <- 2
Q21 <- 2
Q12 <- 2
Q22 <- 2

N <- 1000
library(future.apply)

is_identifiable <- function(pi, rho, alpha) {
    is_col_identifiable <- !any(duplicated(as.vector(alpha %*% t(rho))))
    is_row_identifiable <- !any(duplicated(as.vector(pi %*% alpha)))
    return(is_col_identifiable && is_row_identifiable)
}

plan(multisession(workers = 18))
# With a progressbar from progressr


lbm_param_list <- list()
progressr::with_progress({
    pb <- progressr::progressor(along = 1:N)
    rep_res_lbm <- future_lapply(1:N, function(rep) {
        pi1 <- rdirichlet(1, alpha = rep(1, Q11))

        pi2 <- rdirichlet(1, alpha = rep(1, Q12))

        pi3 <- rdirichlet(1, alpha = rep(1, Q11))

        rho1 <- rdirichlet(1, alpha = rep(1, Q21))

        rho2 <- rdirichlet(1, alpha = rep(1, Q22))

        rho3 <- rdirichlet(1, alpha = rep(1, Q21))

        alpha1 <- matrix(runif(Q11 * Q21, 0.1, 1), nrow = Q11)
        alpha2 <- matrix(runif(Q12 * Q22, 0.1, 1), nrow = Q12)
        alpha3 <- matrix(runif(Q11 * Q21, 0.1, 1), nrow = Q11)
        # Generating graphon distance from all_permutations
        dist_all_perm12 <- colSBM:::graphon_distance_bipartite_all_permutations(
            pis = list(pi1, pi2),
            rhos = list(rho1, rho2),
            alphas = list(alpha1, alpha2)
        )
        dist_all_perm13 <- colSBM:::graphon_distance_bipartite_all_permutations(
            pis = list(pi1, pi3),
            rhos = list(rho1, rho3),
            alphas = list(alpha1, alpha3)
        )
        dist_all_perm23 <- colSBM:::graphon_distance_bipartite_all_permutations(
            pis = list(pi2, pi3),
            rhos = list(rho2, rho3),
            alphas = list(alpha2, alpha3)
        )

        all_perm_dist_vector <- c(dist_all_perm12, dist_all_perm13, dist_all_perm23)
        all_perm_dist_vector <- sqrt(all_perm_dist_vector)

        order_all_perm <- order(all_perm_dist_vector)

        all_perm_tri_ineq <- all(sapply(seq_along(all_perm_dist_vector), function(i) {
            all_perm_dist_vector[i] <= sum(all_perm_dist_vector[-i])
        }))

        # Generating graphon distance from marginals
        dist_marginals12 <- colSBM:::dist_graphon_bipartite_marginals(
            pis = list(pi1, pi2),
            rhos = list(rho1, rho2),
            alphas = list(alpha1, alpha2)
        )
        dist_marginals13 <- colSBM:::dist_graphon_bipartite_marginals(
            pis = list(pi1, pi3),
            rhos = list(rho1, rho3),
            alphas = list(alpha1, alpha3)
        )
        dist_marginals23 <- colSBM:::dist_graphon_bipartite_marginals(
            pis = list(pi2, pi3),
            rhos = list(rho2, rho3),
            alphas = list(alpha2, alpha3)
        )

        marginals_dist_vector <- c(dist_marginals12, dist_marginals13, dist_marginals23)
        marginals_dist_vector <- sqrt(marginals_dist_vector)
        order_marginals <- order(marginals_dist_vector)

        marginal_tri_ineq <- all(sapply(seq_along(marginals_dist_vector), function(i) {
            marginals_dist_vector[i] <= sum(marginals_dist_vector[-i])
        }))


        pb()



        list(df = data.frame(
            rep = rep,
            all_perm_tri_ineq = all_perm_tri_ineq,
            marginal_tri_ineq = marginal_tri_ineq,
            order_kept = all(order_all_perm == order_marginals)
        ), params = list(pi1 = pi1, pi2 = pi2, pi3 = pi3, rho1 = rho1, rho2 = rho2, rho3 = rho3, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3))
    },
    future.seed = TRUE
    )
})

rep_params_lbm <- lapply(rep_res_lbm, function(x) x$params)
rep_res_lbm <- do.call("rbind", lapply(rep_res_lbm, function(x) x$df))

progressr::with_progress({
    pb <- progressr::progressor(along = 1:N)
    rep_res_sbm <- do.call("rbind", future_lapply(1:N, function(rep) {
        pi1 <- rdirichlet(1, alpha = rep(1, Q11))

        pi2 <- rdirichlet(1, alpha = rep(1, Q12))

        pi3 <- rdirichlet(1, alpha = rep(1, Q11))

        alpha1 <- matrix(runif(Q11 * Q21, 0.1, 1), nrow = Q11)
        alpha2 <- matrix(runif(Q12 * Q22, 0.1, 1), nrow = Q12)
        alpha3 <- matrix(runif(Q11 * Q21, 0.1, 1), nrow = Q11)
        # Generating graphon distance from all_permutations
        dist_all_perm12 <- colSBM:::graphon_distance_unipartite_all_permutations(
            pis = list(pi1, pi2),
            alphas = list(alpha1, alpha2)
        )
        dist_all_perm13 <- colSBM:::graphon_distance_unipartite_all_permutations(
            pis = list(pi1, pi3),
            alphas = list(alpha1, alpha3)
        )
        dist_all_perm23 <- colSBM:::graphon_distance_unipartite_all_permutations(
            pis = list(pi2, pi3),
            alphas = list(alpha2, alpha3)
        )

        all_perm_dist_vector <- c(dist_all_perm12, dist_all_perm13, dist_all_perm23)
        all_perm_dist_vector <- sqrt(all_perm_dist_vector)
        order_all_perm <- order(all_perm_dist_vector)

        all_perm_tri_ineq <- sapply(seq_along(all_perm_dist_vector), function(i) {
            all_perm_dist_vector[i] <= sum(all_perm_dist_vector[-i])
        }) |> all()

        # Generating graphon distance from marginals
        dist_marginals12 <- colSBM:::dist_graphon_unipartite_marginals(
            pis = list(pi1, pi2),
            alphas = list(alpha1, alpha2)
        )
        dist_marginals13 <- colSBM:::dist_graphon_unipartite_marginals(
            pis = list(pi1, pi3),
            alphas = list(alpha1, alpha3)
        )
        dist_marginals23 <- colSBM:::dist_graphon_unipartite_marginals(
            pis = list(pi2, pi3),
            alphas = list(alpha2, alpha3)
        )

        marginals_dist_vector <- c(dist_marginals12, dist_marginals13, dist_marginals23)
        marginals_dist_vector <- sqrt(marginals_dist_vector)

        order_marginals <- order(marginals_dist_vector)

        marginal_tri_ineq <- all(sapply(seq_along(marginals_dist_vector), function(i) {
            marginals_dist_vector[i] <= sum(marginals_dist_vector[-i])
        }))

        dist_graphclust12 <- graphclust::sbmNorm(theta1 = list(pi = pi1, gamma = alpha1), theta2 = list(pi = pi2, gamma = alpha2))
        dist_graphclust13 <- graphclust::sbmNorm(theta1 = list(pi = pi1, gamma = alpha1), theta2 = list(pi = pi3, gamma = alpha3))
        dist_graphclust23 <- graphclust::sbmNorm(theta1 = list(pi = pi2, gamma = alpha2), theta2 = list(pi = pi3, gamma = alpha3))

        graphclust_dist_vector <- c(dist_graphclust12, dist_graphclust13, dist_graphclust23)
        graphclust_dist_vector <- sqrt(graphclust_dist_vector)
        order_graphclust <- order(graphclust_dist_vector)

        graphclust_tri_ineq <- all(sapply(seq_along(graphclust_dist_vector), function(i) {
            graphclust_dist_vector[i] <= sum(graphclust_dist_vector[-i])
        }))

        pb()

        data.frame(
            all_perm_tri_ineq = all_perm_tri_ineq,
            marginal_tri_ineq = marginal_tri_ineq,
            graphclust_tri_ineq = graphclust_tri_ineq,
            order_kept = all(order_all_perm == order_marginals),
            order_kept_graphclust = all(order_all_perm == order_graphclust)
        )
    },
    future.seed = TRUE
    ))
})

library(dplyr)
rep_res_sbm |>
    select(all_perm_tri_ineq, marginal_tri_ineq, graphclust_tri_ineq, order_kept, order_kept_graphclust) |>
    summarise_all(mean)

rep_res_lbm |>
    select(all_perm_tri_ineq, marginal_tri_ineq, order_kept) |>
    summarise_all(mean)

non_working_lbm <- rep_res_lbm |>
    filter(!all_perm_tri_ineq) |>
    select(rep) |>
    as.vector() |>
    unlist() |>
    unname()

non_working_lbm_params <- rep_params_lbm[non_working_lbm]

# Verification of nullity for same parameters
progressr::with_progress({
    pb <- progressr::progressor(along = 1:N)
    res_null_lbm <- do.call("rbind", future_lapply(1:N, function(rep) {
        pi1 <- rdirichlet(1, alpha = rep(1, Q11))

        # Shuffle pi1
        order1 <- sample(1:Q11, replace = FALSE)

        rho1 <- rdirichlet(1, alpha = rep(1, Q21))
        order2 <- sample(1:Q21, replace = FALSE)

        alpha1 <- matrix(runif(Q11 * Q21, 0.1, 1), nrow = Q11)

        # Generating graphon distance from all_permutations
        dist_all_perm <- colSBM:::graphon_distance_bipartite_all_permutations(
            pis = list(pi1, pi1[order1]),
            rhos = list(rho1, rho1[order2]),
            alphas = list(alpha1, alpha1[order1, order2])
        )


        # Generating graphon distance from marginals
        dist_marginals <- colSBM:::dist_graphon_bipartite_marginals(
            pis = list(pi1, pi1[order1]),
            rhos = list(rho1, rho1[order2]),
            alphas = list(alpha1, alpha1[order1, order2])
        )


        pb()

        data.frame(
            dist_all_perm = dist_all_perm,
            dist_marginals = dist_marginals,
            marginals_attained_min = dist_marginals == dist_all_perm,
            min_is_zero = all(dist_marginals == 0),
            is_identifiable = is_identifiable(pi1, rho1, alpha1)
        )
    },
    future.seed = TRUE
    ))
})

res_null_lbm |>
    select(dist_all_perm, dist_marginals, marginals_attained_min, min_is_zero) |>
    summarise(
        dist_all_perm = mean(dist_all_perm),
        dist_marginals = mean(dist_marginals),
        marginals_attained_min = sum(marginals_attained_min) / n(),
        min_is_zero = sum(min_is_zero) / n()
    )
