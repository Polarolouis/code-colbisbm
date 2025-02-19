suppressPackageStartupMessages(library("colSBM"))
suppressPackageStartupMessages(library(progressr))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(future.callr))
set.seed(1234)

handlers(global = TRUE)

options(future.globals.maxSize = Inf)

plan(list(
    tweak("callr", workers = floor(parallelly::availableCores(omit = 1L) / 12L)),
    tweak("callr", workers = 4L),
    tweak("callr", workers = 3L)
))

# Network param
nr <- 90
nc <- 90

# Changing parameters
epsilons_pi <- seq(from = 0.0, to = 0.28, by = 0.035)
epsilons_rho <- seq(from = 0.0, to = 0.28, by = 0.035)
pi1 <- matrix(rep(1 / 3, 3), nrow = 1)
rho1 <- matrix(rep(1 / 3, 3), nrow = 1)
ea <- 0.16
alpha <- 0.25 + matrix(
    c(
        3 * ea, 2 * ea, ea,
        2 * ea, 2 * ea, -ea,
        ea,     -ea,    ea
    ),
    byrow = TRUE, nrow = 3, ncol = 3
)


prob_order <- seq(1, 3)
prob_order <- t(sapply(combinat::permn(prob_order), function(v) v))

repetitions <- seq(1, 3)

conditions <- tidyr::crossing(
    epsilon_pi = epsilons_pi,
    epsilon_rho = epsilons_rho,
    pi2_order = prob_order,
    rho2_order = prob_order,
    repetition = repetitions
)

#  Data params
main_dir <- file.path("simulations", "model_selection")

if (!dir.exists(main_dir)) {
    dir.create(main_dir, recursive = TRUE)
}

start_time <- format(Sys.time(), "%d-%m-%Y_%H-%M-%S")
temp_dir <- file.path(main_dir, paste0("tmp", start_time))

if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
}

file_save <- file.path(main_dir, paste0(
    "model_selection", start_time, "_.Rds"
))


tictoc::tic()
row_conditions <- seq_len(nrow(conditions))
results <- future_lapply(row_conditions, function(s) {
    start_time_condition <- Sys.time()
    epsilon_pi <- conditions[s, ]$epsilon_pi
    epsilon_rho <- conditions[s, ]$epsilon_rho

    # Computing the vector with the epsilons
    current_pi2 <- c(
        1 / 3 - epsilon_pi,
        1 / 3,
        1 / 3 + epsilon_pi
    )
    current_rho2 <- c(
        1 / 3 - epsilon_rho,
        1 / 3,
        1 / 3 + epsilon_rho
    )

    current_pi3 <- c(
        1 / 3 + epsilon_pi,
        1 / 3,
        1 / 3 - epsilon_pi
    )
    current_rho3 <- c(
        1 / 3 + epsilon_rho,
        1 / 3,
        1 / 3 - epsilon_rho
    )

    # Permutating the vectors
    current_pi2 <- current_pi2[conditions[s, ]$pi2_order]
    current_rho2 <- current_rho2[conditions[s, ]$rho2_order]
    current_pi3 <- current_pi3[conditions[s, ]$pi2_order]
    current_rho3 <- current_rho3[conditions[s, ]$rho2_order]

    netlist_generated <- list(
        generate_bipartite_collection(
            nr, nc, pi1, rho1,
            alpha,
            M = 1, return_memberships = TRUE
        )[[1]],
        generate_bipartite_collection(
            nr, nc, current_pi2, current_rho2,
            alpha,
            M = 1, return_memberships = TRUE
        )[[1]],
        generate_bipartite_collection(
            nr, nc, current_pi3, current_rho3,
            alpha,
            M = 1, return_memberships = TRUE
        )[[1]]
    )

    # Extracting the incidence matrices
    netlist <- lapply(seq_along(netlist_generated), function(m) {
        return(netlist_generated[[m]]$incidence_matrix)
    })

    # Estimating the models

    list_of_models <- future.apply::future_lapply(c(
        "iid", "pi", "rho", "pirho"
    ), function(model) {
        fitted_bisbmpop <- estimate_colBiSBM(
            netlist = netlist,
            colsbm_model = model,
            nb_run = 3L,
            global_opts = list(
                verbosity = 0,
                plot_details = 0,
                nb_cores = parallelly::availableCores(omit = 1L)
            ),
            fit_opts = list(
                max_vem_steps = 3000L
            )
        )

        return(fitted_bisbmpop)
    }, future.seed = TRUE)

    fitted_bisbmpop_iid <- list_of_models[[1]]

    fitted_bisbmpop_iid$sep_BiSBM$M <- fitted_bisbmpop_iid$M
    sep_BiSBM <- fitted_bisbmpop_iid$sep_BiSBM

    fitted_bisbmpop_pi <- list_of_models[[2]]
    fitted_bisbmpop_rho <- list_of_models[[3]]
    fitted_bisbmpop_pirho <- list_of_models[[4]]

    stop_time_condition <- Sys.time()

    # BICLs
    sep_BICL <- sum(fitted_bisbmpop_iid$sep_BiSBM$BICL)
    iid_BICL <- fitted_bisbmpop_iid$best_fit$BICL
    pi_BICL <- fitted_bisbmpop_pi$best_fit$BICL
    rho_BICL <- fitted_bisbmpop_rho$best_fit$BICL
    pirho_BICL <- fitted_bisbmpop_pirho$best_fit$BICL
    BICLs <- c(sep_BICL, iid_BICL, pi_BICL, rho_BICL, pirho_BICL)

    data_frame_output <- data.frame(
        # The conditions
        epsilon_pi = epsilon_pi,
        epsilon_rho = epsilon_rho,
        pi2 = matrix(current_pi2, nrow = 1),
        rho2 = matrix(current_rho2, nrow = 1),
        pi3 = matrix(current_pi3, nrow = 1),
        rho3 = matrix(current_rho3, nrow = 1),
        repetition = as.numeric(conditions[s, ]$repetition),

        # The results
        ## sep
        sep_BICL = sep_BICL,

        ## iid
        iid_BICL = iid_BICL,
        iid_Q1 = fitted_bisbmpop_iid$best_fit$Q[1],
        iid_Q2 = fitted_bisbmpop_iid$best_fit$Q[2],

        ## pi
        pi_BICL = pi_BICL,
        pi_Q1 = fitted_bisbmpop_pi$best_fit$Q[1],
        pi_Q2 = fitted_bisbmpop_pi$best_fit$Q[2],

        ## pi
        rho_BICL = rho_BICL,
        rho_Q1 = fitted_bisbmpop_rho$best_fit$Q[1],
        rho_Q2 = fitted_bisbmpop_rho$best_fit$Q[2],

        ## pirho
        pirho_BICL = pirho_BICL,
        pirho_Q1 = fitted_bisbmpop_pirho$best_fit$Q[1],
        pirho_Q2 = fitted_bisbmpop_pirho$best_fit$Q[2],

        # Preferred model
        preferred_model = c(
            "sep", "iid",
            "pi", "rho", "pirho"
        )[which.max(BICLs)],
        elapsed_secs = difftime(stop_time_condition,
            start_time_condition,
            units = "sec"
        )
    )
    message("Finished step ", s, "/", nrow(conditions))

    #  Saving temp
    temp_file_save <- file.path(temp_dir, paste0(
        "conditions_", s, "_on_",
        nrow(conditions), ".Rds"
    ))

    saveRDS(object = data_frame_output, file = temp_file_save)

    return(data_frame_output)
}, future.seed = TRUE)
tictoc::toc()
full_data_frame <- do.call(rbind, results)

saveRDS(full_data_frame,
    file = file_save
)
message("Finished simulations.")
