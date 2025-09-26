library(sbm)
library(here)
library(ggplot2) # For plotting
library(patchwork) # For arranging plots
library(dplyr)
library(tikzDevice) # For saving plots as .tex files
library(knitr) # For kable
library(colSBM) # For the colSBM model

# Set up tikzDevice to use standalone document class
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")


incidence_matrices <- readRDS(here("data", "baldock2019-binary-matrices.Rds"))

# Fit the model
set.seed(123)
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

short_names <- c("Bristol", "Edinburgh", "Leeds", "Reading")

if (!dir.exists(here("figures", "applications", "baldock"))) {
    dir.create(here("figures", "applications", "baldock"), recursive = TRUE)
}

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
        theme(axis.ticks = element_blank(), text = element_text(size = 30, family = "serif"))
    return(g)
}

plot_list <- lapply(seq_along(fit_list), function(fit_id) {
    plot_lbm_with_con(fit = fit_list[[fit_id]], title = short_names[fit_id])
})

full_plot <- wrap_plots(plot_list, ncol = 2)


ggsave(here("figures", "applications", "baldock", paste0(
    "baldock2019-fit-sep-BiSBM.pdf")), full_plot)
