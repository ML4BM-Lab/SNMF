
#' Spatially regularized Non-negative Matrix Factorization (sNMF)
#'
#' Performs non-negative matrix factorization (NMF) on count data, incorporating
#' spatial relationships between observations through a spatial similarity matrix `S`.
#'
#' @param counts Numeric matrix of counts (features × observations), e.g., genes × spots.
#' @param S Numeric spatial similarity matrix, e.g., as returned by \code{load_data}.
#' @param k Integer, number of components used for the factorization, e.g. the number of cell-types used for deconvolution.
#' @param Winit Initialization of W (default NULL).
#' @param Hinit Initialization of H (default NULL).
#' @param niter Integer, maximum number of iterations for the NMF algorithm (default 2000).
#' @param Rupdate_iter Integer, number of iterations after which dispersion matrix R is updated (default 10).
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
#'   \item Runs the `SNMF` algorithm on the GPU, with `k` factors.
#'   \item Normalizes the W (basis) and H (coefficient) matrices using the specified quantile.
#'   \item Computes a normalized H matrix (`HC`) representing proportions per spot.
#'   \item Assigns row and column names to W and H matrices based on the input counts.
#' }
#'
#' @return A list with six elements:
#' \describe{
#'   \item{W}{Numeric matrix (features × k), representing basis vectors for features.}
#'   \item{H}{Numeric matrix (k × observations), representing coefficients for each observation.}
#'   \item{phi}{Numeric matrix (features × observations), representing the inverse-dispersion parameter for each entry.}
#'   \item{alpha}{Numeric vector of the number of features as dimension, representing the inverse-dispersion contribution of each gene.}
#'   \item{beta}{Numeric vector of the number of observations as dimension, representing the inverse-dispersion contribution of each spot.}
#'   \item{niter}{Number of iterations run until convergence.}
#' }
#'
#' @examples
#' \dontrun{
#' # counts and S returned from load_data()
#' result <- snmf(counts, S, niter = 1000, Rupdate_iter=20, tol = 1e-4, num_initializations = 5)
#' W <- result$W
#' H <- result$H
#' }
#'
#' @export
snmf <- function(
    counts, 
    S, 
    k, 
    Winit = NULL,
    Hinit = NULL,
    niter=2000, 
    Rupdate_iter=10, 
    tol=1e-4, 
    num_initializations=10, 
    probs=0.75, 
    seed=42
) {

    gpu_counts <-  GPUmatrix::gpu.matrix(counts, dtype = "float32")
    S <-  GPUmatrix::gpu.matrix(S, dtype = "float32")

    set.seed(seed)

    output <- factorize(
        gpu_counts, 
        S = S, 
        k = k,
        Winit = Winit,
        Hinit = Hinit,
        niter = niter, 
        Rupdate_iter = Rupdate_iter,
        tol = tol, 
        num_initializations=num_initializations,
    )

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

    phi <- as.matrix(output$phi)
    rownames(phi) <- rownames(counts)
    colnames(phi) <- colnames(counts)

    alpha <- matrix(output$alpha, nrow = 1)
    beta <- matrix(output$beta, nrow = 1)
    colnames(alpha) <- rownames(counts)
    colnames(beta) <- colnames(counts)

    return(list(
        W=W, 
        H=H, 
        phi=phi, 
        alpha=alpha, 
        beta=beta, 
        niter=output$niter
    ))
}
