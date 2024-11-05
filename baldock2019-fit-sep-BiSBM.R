library(sbm)
library(here)
library(ggplot2) # For plotting
library(patchwork) # For arranging plots
library(tikzDevice) # For saving plots as .tex files

# Set up tikzDevice to use standalone document class
options(tikzDocumentDeclaration = "\\documentclass[10pt]{standalone}")


incidence_matrices <- readRDS(here("data", "baldock2019-binary-matrices.Rds"))

# Fit the model
set.seed(123)
fit_list <- lapply(incidence_matrices, function(incidence_matrix) {
    fit <- estimateBipartiteSBM(
        netMat = incidence_matrix,
        model = "bernoulli", dimLabels = c(row = "pollinators", col = "plants"),
    )
    return(fit)
})

short_names <- c("Bristol", "Edinburgh", "Leeds", "Reading")

if (!dir.exists(here("figures", "applications", "baldock"))) {
    dir.create(here("figures", "applications", "baldock"), recursive = TRUE)
}

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
