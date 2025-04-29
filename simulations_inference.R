necessary_packages <- c(
  "remotes", "tictoc", "combinat", "parallel", "colSBM", "future.apply", "future.callr"
)

options(future.globals.maxSize = Inf)

if (!all((necessary_packages %in% installed.packages()))) {
  install.packages(necessary_packages[-length(necessary_packages)])
}

# Sourcing all necessary files
library(colSBM)
library(future.callr)

future::plan(tweak("callr", workers = parallelly::availableCores(omit = 1L)))
set.seed(1234)

# Network param
nr <- 2 * 120
nc <- 2 * 120
M <- 2

# Changing parameters
base_alpha <- matrix(0.25, nrow = 4, ncol = 4)
epsilon_alpha <- seq(from = 0.0, to = 0.24, by = 0.03)

base_pi1 <- c(0.2, 0.4, 0.4, 0)
rho1 <- rep(0.25, 4)

pi2 <- rep(0.25, 4)
base_rho2 <- c(0, 1 / 3, 1 / 3, 1 / 3)

pi1 <- matrix(unlist(combinat::permn(base_pi1)), byrow = TRUE, ncol = 4L)
pi1 <- pi1[!duplicated(pi1), ]

rho2 <- matrix(unlist(combinat::permn(base_rho2)), byrow = TRUE, ncol = 4L)
rho2 <- rho2[!duplicated(rho2), ]

repetition <- seq.int(3L)

conditions <- tidyr::crossing(epsilon_alpha, pi1, rho2, repetition)

# Filter conditions to prevent the same blocks from being empty
conditions <- conditions[
  !apply(
    conditions[["pi1"]][, 1L:4L] == 0L &
      conditions[["rho2"]][, 1L:4L] == 0L,
    1L, any
  ),
]
#  Data params
main_dir <- file.path("simulations", "inference")

if (!dir.exists(main_dir)) {
  dir.create(main_dir, recursive = TRUE)
}

start_time <- format(Sys.time(), "%d-%m-%Y_%H-%M-%S")
temp_dir <- file.path(main_dir, paste0("tmp", start_time))

if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
}

file_save <- file.path(main_dir, paste0("bernoulli_inference_", start_time, ".Rds"))

missing_conditions_file <- file.path(main_dir, "missing_conditions.Rds")
tictoc::tic()
if (file.exists(missing_conditions_file)) {
  message("Resuming from missing conditions.")
  row_conditions <- readRDS(missing_conditions_file)
} else {
  message("Starting from scratch.")
  row_conditions <- seq_len(nrow(conditions))
}

message("Starting bernoulli inference simulation.")
conditions_rows <- seq_len(nrow(conditions))
results <- future.apply::future_lapply(conditions_rows, function(s) {
  if (!(s %in% row_conditions)) {
    message("Skipping condition ", s, " on ", nrow(conditions))
    return(NULL)
  }
  message(
    "Starting condition ", s, " on ", nrow(conditions),
    " with epsilon_alpha = ", conditions[s, ]$epsilon_alpha,
    " and repetition = ", conditions[s, ]$repetition
  )
  start_time_condition <- Sys.time()
  ea <- conditions[s, ]$epsilon_alpha
  current_pi1 <- conditions[s, ]$pi1
  current_rho2 <- conditions[s, ]$rho2

  current_alpha <- base_alpha + matrix(
    c(
      3 * ea, 2 * ea, ea, -ea,
      2 * ea, 2 * ea, -ea, ea,
      ea, -ea, ea, 2 * ea,
      -ea, ea, 2 * ea, 0
    ),
    byrow = TRUE, nrow = 4, ncol = 4
  )

  # Compute supports
  Cpi1 <- matrix(c(current_pi1, pi2), byrow = TRUE, nrow = 2 * M) > 0
  Cpi2 <- matrix(c(rho1, current_rho2), byrow = TRUE, nrow = 2 * M) > 0

  netlist_generated <- c(
    generate_bipartite_collection(
      nr, nc, conditions[s, ]$pi1, rho1,
      current_alpha,
      M = M, distribution = "bernoulli",
      return_memberships = TRUE
    ),
    generate_bipartite_collection(
      nr, nc, pi2, conditions[s, ]$rho2,
      current_alpha,
      distribution = "bernoulli",
      M = M, return_memberships = TRUE
    )
  )
  netlist <- lapply(seq_along(netlist_generated), function(m) {
    return(netlist_generated[[m]]$incidence_matrix)
  })

  row_clusterings <- lapply(seq_along(netlist_generated), function(m) {
    return(netlist_generated[[m]]$row_blockmemberships)
  })

  col_clusterings <- lapply(seq_along(netlist_generated), function(m) {
    return(netlist_generated[[m]]$col_blockmemberships)
  })

  full_row_clustering <- as.vector(sapply(
    seq.int(2 * M),
    function(m) row_clusterings[[m]]
  ))

  full_col_clustering <- as.vector(sapply(
    seq.int(2 * M),
    function(m) col_clusterings[[m]]
  ))

  fit_opts <- list(max_vem_steps = 5000L)

  fitted_bisbmpop_iid <- estimate_colBiSBM(
    netlist = netlist,
    colsbm_model = "iid",
    nb_run = 3L,
    distribution = "bernoulli",
    global_opts = list(
      verbosity = 0,
      plot_details = 0,
      nb_cores = parallelly::availableCores(omit = 1)
    ),
    fit_opts = fit_opts
  )

  # Handling a problem with sep_BiSBM$M
  fitted_bisbmpop_iid$sep_BiSBM$M <- fitted_bisbmpop_iid$M
  sep_BiSBM <- fitted_bisbmpop_iid$sep_BiSBM

  fitted_bisbmpop_pi <- estimate_colBiSBM(
    netlist = netlist,
    colsbm_model = "pi",
    nb_run = 3L,
    distribution = "bernoulli",
    global_opts = list(
      verbosity = 0,
      plot_details = 0,
      nb_cores = parallelly::availableCores(omit = 1)
    ),
    sep_BiSBM = sep_BiSBM,
    fit_opts = fit_opts
  )

  fitted_bisbmpop_rho <- estimate_colBiSBM(
    netlist = netlist,
    colsbm_model = "rho",
    nb_run = 3L,
    distribution = "bernoulli",
    global_opts = list(
      verbosity = 0,
      plot_details = 0,
      nb_cores = parallelly::availableCores(omit = 1)
    ),
    sep_BiSBM = sep_BiSBM,
    fit_opts = fit_opts
  )

  fitted_bisbmpop_pirho <- estimate_colBiSBM(
    netlist = netlist,
    colsbm_model = "pirho",
    nb_run = 3L,
    distribution = "bernoulli",
    global_opts = list(
      verbosity = 0,
      plot_details = 0,
      nb_cores = parallelly::availableCores(omit = 1)
    ),
    sep_BiSBM = sep_BiSBM,
    fit_opts = fit_opts
  )

  stop_time_condition <- Sys.time()

  ##  Preparing date for export
  # BICLs
  sep_BICL <- sum(sep_BiSBM$BICL)
  iid_BICL <- fitted_bisbmpop_iid$best_fit$BICL
  pi_BICL <- fitted_bisbmpop_pi$best_fit$BICL
  rho_BICL <- fitted_bisbmpop_rho$best_fit$BICL
  pirho_BICL <- fitted_bisbmpop_pirho$best_fit$BICL
  BICLs <- c(sep_BICL, iid_BICL, pi_BICL, rho_BICL, pirho_BICL)

  # ARIs

  compute_mean_ARI <- function(model) {
    # We compute the mean amongst the two networks and return values for
    # rows and columns in a vector
    # sapply ives a matrix with in row the axis ARIs
    # and in columns the networks
    #    1     2
    # ax row1  row2
    # ay col1  col2
    rowMeans(sapply(seq.int(model$M), function(m) {
      matrix(c(
        aricode::ARI(model$Z[[m]][[1]], row_clusterings[[m]]),
        aricode::ARI(model$Z[[m]][[2]], col_clusterings[[m]])
      ), nrow = 2, ncol = 1)
    }))
  }


  compute_double_ARI <- function(model) {
    model_row_Z <- as.vector(sapply(
      seq.int(model$M),
      function(m) model$Z[[m]][[1]]
    ))

    model_col_Z <- as.vector(sapply(
      seq.int(model$M),
      function(m) model$Z[[m]][[2]]
    ))

    return(list(
      aricode::ARI(model_row_Z, full_row_clustering),
      aricode::ARI(model_col_Z, full_col_clustering)
    ))
  }

  sep_mean_ARIs <- compute_mean_ARI(sep_BiSBM)
  iid_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_iid$best_fit)
  pi_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_pi$best_fit)
  rho_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_rho$best_fit)
  pirho_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_pirho$best_fit)

  sep_double_ARIs <- compute_double_ARI(fitted_bisbmpop_iid$sep_BiSBM)
  iid_double_ARIs <- compute_double_ARI(fitted_bisbmpop_iid$best_fit)
  pi_double_ARIs <- compute_double_ARI(fitted_bisbmpop_pi$best_fit)
  rho_double_ARIs <- compute_double_ARI(fitted_bisbmpop_rho$best_fit)
  pirho_double_ARIs <- compute_double_ARI(fitted_bisbmpop_pirho$best_fit)

  data_frame_output <- data.frame(
    # The conditions
    epsilon_alpha = ea,
    pi1 = current_pi1,
    rho2 = current_rho2,
    repetition = as.numeric(conditions[s, 4]),
    # The results
    ## sep
    sep_BICL = sep_BICL,
    sep_mean_row_ARI = sep_mean_ARIs[1],
    sep_mean_col_ARI = sep_mean_ARIs[2],
    sep_double_row_ARI = sep_double_ARIs[[1]],
    sep_double_col_ARI = sep_double_ARIs[[2]],

    ## iid
    iid_BICL = iid_BICL,
    iid_mean_row_ARI = iid_mean_ARIs[1],
    iid_mean_col_ARI = iid_mean_ARIs[2],
    iid_double_row_ARI = iid_double_ARIs[[1]],
    iid_double_col_ARI = iid_double_ARIs[[2]],
    iid_Q1 = fitted_bisbmpop_iid$best_fit$Q[1],
    iid_Q2 = fitted_bisbmpop_iid$best_fit$Q[2],

    ## pi
    pi_BICL = pi_BICL,
    pi_mean_row_ARI = pi_mean_ARIs[1],
    pi_mean_col_ARI = pi_mean_ARIs[2],
    pi_double_row_ARI = pi_double_ARIs[[1]],
    pi_double_col_ARI = pi_double_ARIs[[2]],
    pi_Q1 = fitted_bisbmpop_pi$best_fit$Q[1],
    pi_Q2 = fitted_bisbmpop_pi$best_fit$Q[2],

    ## pi
    rho_BICL = rho_BICL,
    rho_mean_row_ARI = rho_mean_ARIs[1],
    rho_mean_col_ARI = rho_mean_ARIs[2],
    rho_double_row_ARI = rho_double_ARIs[[1]],
    rho_double_col_ARI = rho_double_ARIs[[2]],
    rho_Q1 = fitted_bisbmpop_rho$best_fit$Q[1],
    rho_Q2 = fitted_bisbmpop_rho$best_fit$Q[2],

    ## pirho
    pirho_BICL = pirho_BICL,
    pirho_mean_row_ARI = pirho_mean_ARIs[1],
    pirho_mean_col_ARI = pirho_mean_ARIs[2],
    pirho_double_row_ARI = pirho_double_ARIs[[1]],
    pirho_double_col_ARI = pirho_double_ARIs[[2]],
    pirho_Q1 = fitted_bisbmpop_pirho$best_fit$Q[1],
    pirho_Q2 = fitted_bisbmpop_pirho$best_fit$Q[2],
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

  #  Saving inhabitual data
  if (all(unlist(pirho_mean_ARIs) == 1L) & any(unlist(pirho_double_ARIs) < 1L)) {
    warning("Incorrect result encountered, saving.")
    incorrect_filepath <- file.path(temp_dir, paste0(
      "incorrect_conditions_", s, "_on_",
      nrow(conditions), ".Rds"
    ))

    inc_data <- list(
      epsilon_alpha = ea,
      pi1 = current_pi1,
      rho2 = current_rho2,
      alpha = current_alpha,
      pirho_double_row_ARI = pirho_double_ARIs[[1]],
      pirho_double_col_ARI = pirho_double_ARIs[[2]],
      netlist = netlist_generated
    )
    saveRDS(object = inc_data, file = incorrect_filepath)
  }

  return(data_frame_output)
},
future.seed = TRUE
)


tictoc::toc()
full_data_frame <- do.call(rbind, results)

saveRDS(full_data_frame,
  file = file_save
)
message("Finished simulations.")
