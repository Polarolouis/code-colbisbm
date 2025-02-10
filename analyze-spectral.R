devtools::load_all("../colSBM/")

# Load the data
data("dorebipartite")

set.seed(123L)
mesure_sim <- function(nstart, iter.max, rep.kmeans = 100) {
    # Compare the spectral clustering of the two models
    mat <- dorebipartite[[1]] %*% t(dorebipartite[[1]])
    clust_list <- lapply(1:rep.kmeans, function(i) {
        spectral_clustering(X = mat, K = 12, kmeans.nstart = nstart, kmeans.iter.max = iter.max)
    })

    ari_dist <- function(clust1, clust2) {
        ari <- aricode::ARI(clust1, clust2)
        return(ari)
    }

    ari_mat <- outer(clust_list, clust_list, Vectorize(ari_dist))

    mean(ari_mat[lower.tri(ari_mat)])
}

mat <- dorebipartite[[1]] %*% t(dorebipartite[[1]])

km <- microbenchmark(
    "n10" = {
        kmeans(mat, centers = 20, iter.max = 50, nstart = 10)
    },
    "n100" = {
        kmeans(mat, centers = 20, iter.max = 50, nstart = 100)
    },
    "n400" = {
        kmeans(mat, centers = 20, iter.max = 50, nstart = 400)
    }
)
