library("rlang")

check_identifiability_iid <- function(nr, nc, Q, alpha, pi, rho, call = caller_env()) {
    check_required(nr)
    check_required(nc)
    check_required(Q)
    check_required(alpha)
    check_required(pi)
    check_required(rho)

    if (all(nr < (2 * Q[2] - 1))) {
        cli::cli_abort(c("x" = "No networks in the collection has a sufficient number of row nodes (n1) for Q2 = {.val {Q[2]}}", "i" = "One of the networks must have n1 >= {.val {2*Q[2]-1}}, current max is {.val {max(nr)}}"))
    }
    if (all(nc < (2 * Q[1] - 1))) {
        cli::cli_abort(c("x" = "No networks in the collection has a sufficient number of column nodes for (n2) Q1 = {.val {Q[1]}}", "i" = "One of the networks must have n2 >= {.val {2*Q[1]-1}}, current max is {.val {max(nc)}}"))
    }

    if (any(duplicated(as.vector(alpha %*% rho)))) {
        cli::cli_abort(c("x" = "All values of {.code alpha %*% rho} must be unique !"))
    }

    if (any(duplicated(as.vector(t(pi) %*% alpha)))) {
        cli::cli_abort(c("x" = "All values of {.code t(pi) %*% alpha} must be unique !"))
    }

    cli::cli_alert_success("Identifiability for iid model passed !")
}

check_identifiability_pirho <- function(nr, nc, Q, alpha, pis, rhos, call = caller_env()) {
    check_required(nr)
    check_required(nc)
    check_required(Q)
    check_required(alpha)
    check_required(pis)
    check_required(rhos)

    M <- length(nr)

    if (any(nr < (2 * Q[2] - 1))) {
        cli::cli_abort(c("x" = "All networks in the collection must have a sufficient number of row nodes (n1) for Q2 = {.val {Q[2]}}", "i" = "They must have n1 >= {.val {2*Q[2]-1}}, current min is {.val {min(nr)}}"))
    }
    if (any(nc < (2 * Q[1] - 1))) {
        cli::cli_abort(c("x" = "All networks in the collection must have a sufficient number of column nodes for (n2) Q1 = {.val {Q[1]}}", "i" = "They must have n2 >= {.val {2*Q[1]-1}}, current min is {.val {min(nc)}}"))
    }

    alpha_ms <- lapply(seq(M), function(m) {
        alpha[pis[[m]] > 0, rhos[[m]] > 0]
    })

    row_identif_vec <- sapply(seq(M), function(m) {
        any(duplicated(as.vector(alpha_ms[[m]] %*% rhos[[m]])))
    })
    if (any(row_identif_vec)) {
        net_fail_row_identif <- which(row_identif_vec)
        cli::cli_abort(c(
            "x" = "All marginals on {.code alpha[[m]] %*% rho[[m]]} must be different.",
            "i" = "Currently {cli::no(cli::qty(length(net_fail_row_identif)))} network{?s} do not respect the condition: {.val {net_fail_row_identif}}."
        ))
    }

    col_identif_vec <- sapply(seq(M), function(m) {
        any(duplicated(as.vector(t(pis[[m]]) %*% alpha_ms[[m]])))
    })
    if (any(col_identif_vec)) {
        net_fail_col_identif <- which(col_identif_vec)
        cli::cli_abort(c(
            "x" = "All marginals on {.code t(pis[[m]]) %*% alpha[[m]]} must be different.",
            "i" = "Currently {cli::no(cli::qty(length(net_fail_col_identif)))} network{?s} do not respect the condition: {.val {net_fail_col_identif}}."
        ))
    }
    if (any(duplicated(as.vector(alpha)))) {
        cli::cli_abort(c("x" = "All entries of `alpha` must be unique !"))
    }
    cli::cli_alert_success("Identifiability for pirho model passed !")
}
