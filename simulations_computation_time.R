library("colSBM")
library("future")
library("futurize")
library("future.callr")
source("utils.R")

save_path <- "simulations/computation_time"

if (!dir.exists(save_path)) {
    dir.create(save_path)
}
epoch <- as.integer(Sys.time())

args <- commandArgs(trailingOnly = TRUE)

model <- "iid"
TEST_NULL_PROP <- FALSE

if (length(args) == 2) {
    model <- args[1]
    TEST_NULL_PROP <- ifelse(args[2], TRUE, FALSE)
} else {
    message("No or incorrect sim params given, defaulting to model : ", model, " and", ifelse(TEST_NULL_PROP, "", " not"), " testing with null block props.")
}

nb_rep <- 3
max_Q <- 8
M <- 3
nr <- 100
nc <- 100

base_filename <- paste0("computation_time_model_", model, "_Qmax_", max_Q, "_M_", M, "_nr_", nr, "_nc_", nc, "_", epoch)

save_file <- file.path(save_path, paste0(base_filename, ".Rds"))
temp_path <- file.path(save_path, paste0("tmp", base_filename))
if (!dir.exists(temp_path)) {
    dir.create(temp_path)
}

seq_Q <- seq(2, max_Q)

build_mod_alpha <- function(Q) {
    alpha <- matrix(0, nrow = Q[1], ncol = Q[2])
    seq_Q1 <- seq(Q[1])
    seq_Q2 <- seq(Q[2])
    diag(alpha) <- exp(-1 / ifelse(Q[1] <= Q[2], seq_Q1, seq_Q2))

    alpha <- alpha + abs(rnorm(n = Q[1] * Q[2], sd = 0.001))

    return(alpha)
}

conditions <- expand.grid(
    Q1 = seq_Q,
    Q2 = seq_Q,
    rep = seq(nb_rep)
)
set.seed(345)
datasets_iid <- lapply(seq_len(nrow(conditions)), function(row_idx) {
    Q1 <- conditions[row_idx, "Q1"]
    Q2 <- conditions[row_idx, "Q2"]
    rep <- conditions[row_idx, "rep"]
    message(c("Generating dataset for Q = (", Q1, ",", Q2, "), rep = ", rep))

    pis <- exp(-1 / seq(Q1))
    pis <- pis / sum(pis)
    rhos <- exp(-1 / rev(seq(Q2)))
    rhos <- rhos / sum(rhos)
    Q <- c(Q1, Q2)
    alpha_iid <- build_mod_alpha(Q)

    check_identifiability_iid(nr = nr, nc = nc, Q = c(Q1, Q2), alpha = alpha_iid, pi = pis, rho = rhos)

    generate_bipartite_collection(nr = nr, nc = nc, pi = pis, rho = rhos, alpha = alpha_iid, M = M)
})

datasets_pirho <- lapply(seq_len(nrow(conditions)), function(row_idx) {
    Q1 <- conditions[row_idx, "Q1"]
    Q2 <- conditions[row_idx, "Q2"]
    rep <- conditions[row_idx, "rep"]
    message(c("Generating dataset for Q = (", Q1, ",", Q2, "), rep = ", rep))

    pi <- exp(-1 / seq(Q1))
    if (TEST_NULL_PROP && Q1 > 2) {
        # If more than 2 clusters introducing a null block
        pi[length(pi)] <- 0
    }
    pi <- pi / sum(pi)
    rho <- exp(-1 / rev(seq(Q2)))
    if (TEST_NULL_PROP && Q2 > 2) {
        # If more than 2 clusters introducing a null block
        rho[1] <- 0
    }
    rho <- rho / sum(rho)
    Q <- c(Q1, Q2)
    alpha_pirho <- build_mod_alpha(Q)


    pis <- lapply(seq(M), function(m) sample(pi))
    rhos <- lapply(seq(M), function(m) sample(rho))

    collection <- generate_bipartite_collection(nr = nr, nc = nc, pi = pis, rho = rhos, alpha = alpha_pirho, M = M, model = "pirho")


    check_identifiability_pirho(nr = nr, nc = nc, Q = c(Q1, Q2), alpha = alpha_pirho, pis = pis, rhos = rhos)

    collection
})

datasets <- datasets_iid

if (model != "iid") {
    message("Using pirho datasets")
    datasets <- datasets_pirho
} else {
    message("Using iid datasets")
}
options("future.globals.maxSize" = Inf)
nbCores <- parallelly::availableCores(omit = 1L)
nb_run <- 3L
plan(tweak("multisession", workers = nbCores))

results_list <- lapply(seq_len(nrow(conditions)), function(row_idx) {
    Q1 <- conditions[row_idx, "Q1"]
    Q2 <- conditions[row_idx, "Q2"]
    rep <- conditions[row_idx, "rep"]
    message(c("Starting sim for Q = (", Q1, ",", Q2, "), rep = ", rep, ", model : ", model))

    collection <- datasets[[row_idx]]
    start_time <- Sys.time()
    fit_res <- estimate_colBiSBM(netlist = collection, colsbm_model = model, global_opts = list(verbosity = 3L, backend = "no_mc"), nb_run = nb_run)
    time_taken <- difftime(Sys.time(), start_time, units = "secs")

    res_df <- data.frame(model = model, Q1 = Q1, Q2 = Q2, rep = rep, duration = time_taken, correct_Q1 = fit_res$best_fit$Q[1] == Q1, correct_Q2 = fit_res$best_fit$Q[2] == Q2)

    saveRDS(res_df, file = file.path(temp_path, paste0("condition_", row_idx, "_on_", nrow(conditions), ".Rds")))

    return(res_df)
}) |> futurize(seed = TRUE)

results_df <- do.call("rbind", results_list)

saveRDS(results_df, save_file)
