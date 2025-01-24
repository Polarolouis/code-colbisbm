library(colSBM)
library(future.apply)
library(future.callr)
library(here)
library(progressr)
# Load baldock2019 data
incidence_matrices <- readRDS(here("data", "baldock2019-binary-matrices.Rds"))


if (!dir.exists(here("figures", "applications", "baldock"))) {
    dir.create(here("figures", "applications", "baldock"), recursive = TRUE)
}

if (!file.exists(here("data", "baldock-iid-pi-rho-pirho.Rds"))) {
    options(future.globals.maxSize = Inf)
    set.seed(123L)
    plan(list(tweak("callr", workers = 4L), tweak("callr", workers = 3L)))
    # Fit the 4 models : iid, pi, rho, pi-rho in parallel
    models <- c("iid", "pi", "rho", "pirho")
    with_progress({
        p <- progressor(along = models)
        fit_list <- future_lapply(models, function(model) {
            fit <- estimate_colBiSBM(
                netlist = incidence_matrices,
                colsbm_model = model,
                net_id = names(incidence_matrices),
                distribution = "bernoulli",
                fit_opts = list(
                    algo_ve = "fp",
                    minibatch = TRUE,
                    verbosity = 0
                ),
                global_opts = list(
                    verbosity = 0,
                    backend = "future"
                ),
                nb_run = 3L
            )
            p(sprintf("Fitted %s model", model))
            return(fit)
        }, future.seed = TRUE)
    })
    saveRDS(fit_list, here("data", "baldock-iid-pi-rho-pirho.Rds"))
} else {
    fit_list <- readRDS(here("data", "baldock-iid-pi-rho-pirho.Rds"))
}
