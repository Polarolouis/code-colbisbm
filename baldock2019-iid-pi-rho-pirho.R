library(colSBM)
library(future.apply)
library(future.callr)
library(here)
library(progressr)
# library(animint2)
handlers(global = TRUE)
handlers("cli")
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

fit_list <- lapply(fit_list, function(fit) {
    fit$best_fit$net_id <- c("Bristol", "Edinburgh", "Leeds", "Reading")
    fit
})

library(ggplot2)

models <- c("iid", "pi", "rho", "pirho")

lapply(seq_along(models), function(id) {
    pdf(file = here("figures", "applications", "baldock", paste0("baldock2019-joint-", models[id], ".pdf")))
    print(plot(fit_list[[id]]$best_fit, type = "meso", mixture = TRUE, values = TRUE) +
        ggtitle(models[id])) + theme_minimal()
    dev.off()
})

patchwork::wrap_plots(lapply(c("BICL", "ICL", "vbound"), function(criterion) {
    plot(fit_list[[1]], criterion = criterion) + ggtitle(criterion)
}), ncol = 2)
