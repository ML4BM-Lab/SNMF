
#' Spatially regularized Non-negative Matrix Factorization (sNMF)
#'
#' Performs non-negative matrix factorization (NMF) on count data, incorporating
#' spatial relationships between observations through a spatial similarity matrix `S`.
#'
#' @param counts Numeric matrix of counts (features × observations), e.g., genes × spots.
#' @param S Numeric spatial similarity matrix, e.g., as returned by \code{load_data}.
#' @param k Integer, number of components used for the factorization, e.g. the number of cell-types used for deconvolution.
#' @param niter Integer, maximum number of iterations for the NMF algorithm (default 2000).
#' @param tol Numeric, convergence tolerance for the NMF algorithm (default 1e-4).
#' @param num_initializations Integer, number of random initializations for NMF (default 10).
#' @param probs Numeric between 0 and 1, quantile used for normalizing W and H matrices (default 0.75).
#' @param seed Integer, random seed for reproducibility (default 42).
#'
#' @details
#' This function performs spatially regularized NMF on the input count matrix:
#' \enumerate{
#'   \item Converts `counts` and `S` to GPU matrices for efficient computation.
#'   \item Sets the random seed for reproducibility.
#'   \item Runs the `NMFKLMixing` algorithm on the GPU, with `k` factors.
#'   \item Normalizes the W (basis) and H (coefficient) matrices using the specified quantile.
#'   \item Computes a normalized H matrix (`HC`) representing proportions per spot.
#'   \item Assigns row and column names to W and H matrices based on the input counts.
#' }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{W}{Numeric matrix (features × k), representing basis vectors for features.}
#'   \item{H}{Numeric matrix (k × observations), representing coefficients for each observation.}
#' }
#'
#' @examples
#' \dontrun{
#' # counts and S returned from load_data()
#' result <- snmf(counts, S, niter = 1000, tol = 1e-4, num_initializations = 5)
#' W <- result$W
#' H <- result$H
#' }
#'
#' @export
snmf <- function(counts, S, k, niter=2000, tol=1e-4, num_initializations=10, probs=0.75, seed=42) {

    gpu_counts <- gpu.matrix(counts, dtype = "float32")
    S <- gpu.matrix(S, dtype = "float32")

    set.seed(seed)

    output <- NMFKLMixing(gpu_counts, S = S, k = k,
                        niter = niter, tol = tol, num_initializations=num_initializations)

    W <- as.matrix(output$W)
    H <- as.matrix(output$H)

    D <- diag(matrixStats::colQuantiles(W, probs = probs, na.rm = T))
    D_1 <- diag(1/matrixStats::colQuantiles(W, probs = probs, na.rm = T))

    # Normalize the W and H matrices
    W <- W %*% D_1; 
    H <- D %*% H
    HC <- as.matrix(H %*% S)

    # To get an H matrix of proportions.
    HC <- t(t(HC)/colSums(HC))

    # Name spots and genes
    colnames(HC) <- colnames(counts)
    HC <- t(HC)

    rownames(W) <- rownames(counts)

    return(list(W=W, H=H))
}
