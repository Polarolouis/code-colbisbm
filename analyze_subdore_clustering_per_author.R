library(here)
library(ggplot2)
library(patchwork)

data_path <- here("applications", "subdore", "tmp20250418103624")

author_vector <- c("Baldock", "Traveset", "Souza", "Cordeniz", "Trojelsgaard", "Gibson")
model_vector <- c("iid", "pi", "rho", "pirho")


flist <- list.files(data_path)
author_order <- order(factor(sub(".*_(Baldock|Traveset|Souza|Cordeniz|Trojelsgaard|Gibson)_.*", "\\1", flist), levels = author_vector))
flist[author_order]

# Create named nested list
# $Author$iid
# $Author$pi
# $Author$rho
# $Author$pirho

setNames(lapply(author_vector, function(author) {
    # Filter flist by author
    author_files <- flist[grepl(author, flist)]
    # Filte by model
    setNames(lapply(model_vector, function(model) {
        # Filter author_files by model
        model_files <- author_files[grepl(model, author_files)]
        # Read files
        lapply(model_files, function(file) {
            readRDS(file.path(data_path, file))
        })
    }), paste0(model_vector, seq(1, 4)))
}), author_vector) -> nested_list
