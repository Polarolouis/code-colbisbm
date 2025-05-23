# devtools::load_all("../colSBM")
library(colSBM)
library(here)
library(stringr)
library(future.callr)
library(tidyverse)

options(future.globals.maxSize = 8912896000) # 8.3 GB
# matrices <- readRDS(here("data", "dore-binary-matrices.Rds"))
# baldock_matrices <- matrices[grepl("Baldock", x = names(matrices))]
# souza_matrices <- matrices[grepl("Souza", x = names(matrices))]

eps <- 0.3
alpha_assortative <- matrix(0.3, nrow = 3, ncol = 3) +
    matrix(
        c(
            eps, -0.5 * eps, -0.5 * eps,
            -0.5 * eps, eps, -0.5 * eps,
            -0.5 * eps, -0.5 * eps, eps
        ),
        nrow = 3, byrow = TRUE
    )

alpha_core_periphery <- matrix(0.3, nrow = 3, ncol = 3) +
    matrix(
        c(
            1.5 * eps, eps, 0.5 * eps,
            eps, 0.5 * eps, 0,
            0.5 * eps, 0, -0.5 * eps
        ),
        nrow = 3, byrow = TRUE
    )

alpha_disassortative <- matrix(0.3, nrow = 3, ncol = 3) +
    matrix(
        c(
            -0.5 * eps, eps, eps,
            eps, -0.5 * eps, eps,
            eps, eps, -0.5 * eps
        ),
        nrow = 3, byrow = TRUE
    )

nr <- 100
nc <- 100

pi <- c(0.3, 0.5, 0.2)
rho <- c(0.2, 0.5, 0.3)
M <- 3
{
    set.seed(123)
    matrices <- c(
        generate_bipartite_collection(
            nr = nr,
            nc = nc,
            pi = pi,
            rho = rho,
            alpha = alpha_assortative,
            M = M,
            model = "iid"
        ),
        generate_bipartite_collection(
            nr = nr,
            nc = nc,
            pi = pi,
            rho = rho,
            alpha = alpha_core_periphery,
            M = M,
            model = "iid"
        ),
        generate_bipartite_collection(
            nr = nr,
            nc = nc,
            pi = pi,
            rho = rho,
            alpha = alpha_disassortative,
            M = M,
            model = "iid"
        )
    )
}

names(matrices) <- paste0(rep(c("assortative", "core_periphery", "disassortative"), each = 3), ".", seq(1, 3))
if (!file.exists(here("investigate", "iid-penalty0.Rds"))) {
    if (!dir.exists(here("investigate"))) {
        dir.create(here("investigate"), recursive = TRUE)
    }
    plan("callr", workers = 3L)
    fit_no_penalty <- estimate_colBiSBM(
        netlist = matrices,
        colsbm_model = "iid",
        net_id = names(matrices),
        fit_opts = list(penalty_factor = 0),
        global_opts = list(
            verbosity = 0,
            backend = "future"
        )
    )
    saveRDS(fit_no_penalty, here("investigate", "iid-penalty0.Rds"))
} else {
    fit_no_penalty <- readRDS(here("investigate", "iid-penalty0.Rds"))
}
pdf("tutu.pdf")
for (q1 in seq(1, fit_no_penalty$global_opts$Q1_max)) {
    for (q2 in seq(1, fit_no_penalty$global_opts$Q2_max)) {
        model <- fit_no_penalty$model_list[[q1, q2]]
        if (!is.null(model)) {
            image(colSBM:::compute_dissimilarity_matrix(model), main = toString(c(q1, q2)))
        }
    }
}
dev.off()


clusts <- list()
clusts_df <- data.frame()
for (q1 in seq(1, fit_no_penalty$global_opts$Q1_max)) {
    for (q2 in seq(1, fit_no_penalty$global_opts$Q2_max)) {
        model <- fit_no_penalty$model_list[[q1, q2]]
        if (!is.null(model)) {
            clustering <- cutree(hclust(as.dist(colSBM:::compute_dissimilarity_matrix(model)), method = "ward.D2"), k = 2)
            clusts[[paste0(toString(c(q1, q2)))]] <- clustering
            clusts_df <- rbind(clusts_df, as.data.frame(list("clustering" = t(clustering))))
        }
    }
}

outer(clusts, clusts, FUN = Vectorize(function(x, y) {
    aricode::ARI(x, y)
})) -> ari_models

ari_models %>%
    reshape2::melt() -> my_df

library(ggokabeito)
ggplot(my_df) +
    geom_tile(aes(x = Var1, y = Var2, fill = as.factor(value)), show.legend = T) +
    scale_fill_okabe_ito()

unique_clusterings <- unique(clusts)
plan("callr", workers = min(parallelly::availableCores(omit = 1L), length(unique_clusterings)))
lapply(unique_clusterings, function(clustering) {
    lapply(c(1, 2), function(k) {
        fit <- estimate_colBiSBM(
            netlist = matrices[clustering == k],
            colsbm_model = "iid",
            net_id = names(matrices)[clustering == k],
            global_opts = list(
                verbosity = 0,
                backend = "no_mc"
            ),
        )
    })
})
