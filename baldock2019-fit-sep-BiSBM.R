library(sbm)
library(here)
library(ggplot2) # For plotting
library(patchwork) # For arranging plots
<<<<<<< HEAD
library(dplyr)
library(tikzDevice) # For saving plots as .tex files
library(knitr) # For kable
library(colSBM) # For the colSBM model
=======
library(tikzDevice) # For saving plots as .tex files
>>>>>>> Adding first data

# Set up tikzDevice to use standalone document class
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")


incidence_matrices <- readRDS(here("data", "baldock2019-binary-matrices.Rds"))

# Fit the model
set.seed(123)
<<<<<<< HEAD
fit_list <- lapply(seq_along(incidence_matrices), function(idx) {
    # fit <- estimateBipartiteSBM(
    #     netMat = incidence_matrix,
    #     model = "bernoulli", dimLabels = c(row = "pollinators", col = "plants"),
    #     estimOptions = list(
    #         verbosity = 0,
    #         plot = FALSE
    #     )
    # )
    incidence_matrix <- incidence_matrices[[idx]]
    net_id <- names(incidence_matrices)[idx]
    fit <- estimateBipartiteSBM(netMat = incidence_matrix, model = "bernoulli", dimLabels = c(row = "pollinators", col = "plants"), estimOptions = list(verbosity = 0, plot = FALSE))
    refit <- colSBM:::fitBipartiteSBMPop$new(
        A = list(incidence_matrix),
        Q = fit$nbBlocks,
        Z = list(fit$memberships),
        net_id = net_id,
        init_method = "given",
        distribution = "bernoulli",
        free_mixture_row = FALSE,
        free_mixture_col = FALSE,
        fit_opts = list(
            algo_ve = "fp",
            minibatch = TRUE, verbosity = 0
        )
    )
    refit$optimize()
    return(refit)
=======
fit_list <- lapply(incidence_matrices, function(incidence_matrix) {
    fit <- estimateBipartiteSBM(
        netMat = incidence_matrix,
        model = "bernoulli", dimLabels = c(row = "pollinators", col = "plants"),
    )
    return(fit)
>>>>>>> Adding first data
})

short_names <- c("Bristol", "Edinburgh", "Leeds", "Reading")

if (!dir.exists(here("figures", "applications", "baldock"))) {
    dir.create(here("figures", "applications", "baldock"), recursive = TRUE)
}

<<<<<<< HEAD
plot_lbm_with_con <- function(fit, m = 1, title = NULL, oRow = NULL, oCol = NULL) {
    if (is.null(oRow)) {
        oRow <- order(fit$alpha %*% matrixStats::rowMeans2(sapply(fit$pi, function(pi) pi[[2]])), decreasing = TRUE)
    }
    if (is.null(oCol)) {
        oCol <- order(matrixStats::rowMeans2(sapply(fit$pi, function(pi) pi[[1]])) %*% fit$alpha, decreasing = TRUE)
    }
    # Here we extract the adjacency matrix and putting it in a geom_tile compatible format
    plot_df <- as.matrix(fit$A[[m]])[order(oRow[fit$Z[[m]]$row]), order(oCol[fit$Z[[m]]$col])] |>
        reshape2::melt() |>
        # Below we add the alpha connectivity values
        mutate(con = fit$alpha[fit$Z[[m]]$row, fit$Z[[m]]$col][order(oRow[fit$Z[[m]]$row]), order(oCol[fit$Z[[m]]$col])] |>
            reshape2::melt() |>
            select(value) |>
            purrr::as_vector())
    # Here we plot the adjacency matrix
    g <- ggplot(plot_df, aes(x = Var2, y = Var1, fill = value, alpha = value)) +
        geom_tile(aes(x = Var2, y = Var1, alpha = con), fill = "red", linewidth = 0, show.legend = FALSE) +
        geom_tile(show.legend = FALSE) +
        geom_hline(yintercept = cumsum(table(fit$Z[[m]]$row)[rev(oRow)]) + .5, col = "red", size = .3) +
        geom_vline(xintercept = cumsum(table(fit$Z[[m]]$col)[oCol]) + .5, col = "red", size = .3) +
        scale_fill_gradient(low = "white", high = "black") +
        ylab("Pollinators") +
        xlab("Plants") +
        scale_alpha(range = c(0, 1)) +
        scale_x_discrete(breaks = "") +
        scale_y_discrete(breaks = "", guide = guide_axis(angle = 0), limits = rev) +
        coord_equal(expand = FALSE) +
        theme_bw(base_size = 20, base_rect_size = 1, base_line_size = 1) +
        ggtitle(title) +
        theme(axis.ticks = element_blank(), text = element_text(size = 30))
    return(g)
}


=======
>>>>>>> Adding first data
# Save the plots
lapply(seq_along(fit_list), function(i) {
    pdf(
        here("figures", "applications", "baldock", paste0(
            "baldock2019-fit-sep-BiSBM-",
            short_names[i],
            ".pdf"
        )),
        family = "Times"
    )
<<<<<<< HEAD

    print(plot_lbm_with_con(fit_list[[i]], title = short_names[i]))
    dev.off()
})


# Retrieving clusters sizes per network
clustersizes <- sapply(fit_list, function(fit) {
    return(fit$nbBlocks)
})
# Rename rows
# row.names(clustersizes) <- c("$\\hat{Q_1^m}$ Pollinators", "$\\hat{Q_2^m}$ Plants")
# Prepare a kable for latex inclusion
# kable(clustersizes,
#     caption = "Results for separate BiSBM on \\cite{baldockSystemsApproachReveals2019}",
#     col.names = short_names, format = "latex",
#     escape = FALSE
# ) |>
#     cat()

# Graphon distance matrix
dist_mat <- outer(fit_list, fit_list, Vectorize(function(fit1, fit2) {
    param1 <- fit1$parameters
    param2 <- fit2$parameters

    colSBM:::dist_graphon_symmetrization(
        pis = list(param1$pi[[1]], param2$pi[[1]]),
        rhos = list(param1$rho[[1]], param2$rho[[1]]),
        alphas = list(param1$alpha, param2$alpha)
    )
}))

dist_mat_marg <- outer(fit_list, fit_list, Vectorize(function(fit1, fit2) {
    param1 <- fit1$parameters
    param2 <- fit2$parameters

    colSBM:::dist_graphon_marginals(
        pis = list(param1$pi[[1]], param2$pi[[1]]),
        rhos = list(param1$rho[[1]], param2$rho[[1]]),
        alphas = list(param1$alpha, param2$alpha)
    )
}))

# For version
dist_mat_perm <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
    for (j in 1:i) {
        fit1 <- fit_list[[i]]
        fit2 <- fit_list[[j]]
        param1 <- fit1$parameters
        param2 <- fit2$parameters
        dist_mat_perm[i, j] <- colSBM:::dist_graphon_marginals(
            pis = list(param1$pi[[1]], param2$pi[[1]]),
            rhos = list(param1$rho[[1]], param2$rho[[1]]),
            alphas = list(param1$alpha, param2$alpha)
        )
    }
}

dist_mat_marginals <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
    for (j in 1:i) {
        fit1 <- fit_list[[i]]
        fit2 <- fit_list[[j]]
        param1 <- fit1$parameters
        param2 <- fit2$parameters
        dist_mat_marginals[i, j] <- colSBM:::graphon_distance_identif(
            pis = list(param1$pi[[1]], param2$pi[[1]]),
            rhos = list(param1$rho[[1]], param2$rho[[1]]),
            alphas = list(param1$alpha, param2$alpha)
        )
    }
}
=======
    print(plot(fit_list[[i]], type = "data") +
        labs(x = "Plants", y = "Pollinators") +
        ggtitle(short_names[i]) +
        theme(
            strip.text.y = element_blank(),
            strip.text.x = element_blank(),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1
            ),
            text = element_text(size = 30)
        ))
    dev.off()
})
>>>>>>> Adding first data
