library("here")
library("sbm")

# Load the data
subject_matrices_path <- here("dataset/NITRC-ADHD200/brains")

# Read each edge list file and convert it back to a matrix
netlist <- lapply(list.files(subject_matrices_path, pattern = "edge_list", full.names = TRUE), function(file) {
    edge_list <- read.csv(file, header = FALSE)
    adj_mat <- as.matrix(table(edge_list$V1, edge_list$V2))
    oRow <- order(as.numeric(head(rownames(adj_mat), -1)))
    oCol <- order(as.numeric(head(colnames(adj_mat), -1)))
    # Reorder the matrix
    adj_mat <- adj_mat[oRow, oCol]
    class(adj_mat) <- "matrix"
    adj_mat
})
names(netlist) <- gsub(list.files(subject_matrices_path, pattern = "edge_list", full.names = FALSE), pattern = "edge_list_", replacement = "") |> gsub(pattern = ".csv", replacement = "")

fit_sbm <- estimateBipartiteSBM(netMat = adj_mat)
