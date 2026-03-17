
#' Load and preprocess count data for NMF
#'
#' This function filters and preprocesses a count matrix for non-negative
#' matrix factorization (NMF), and computes a spatial similarity matrix `S`
#' based on the positions of the columns.
#'
#' @param counts Numeric matrix of counts, where rows are features (e.g., genes)
#'   and columns are observations (e.g., spots or cells). Column names should
#'   encode positions in the form "x<coord>x<coord>".
#' @param filter_th Numeric, minimum sum per row to keep a feature (default 10).
#'   Features with row sums ≤ `filter_th` are removed.
#' @param tau Numeric, target mean value for optimization of gamma (default 0.5).
#'   This is used in computing the spatial similarity matrix.
#' @param S_th Numeric, threshold below which entries of the similarity matrix
#'   are set to zero (default 1e-3).
#'
#' @details
#' The function performs the following steps:
#' 1. Filters rows of the counts matrix with low total counts.
#' 2. Extracts x and y positions from column names.
#' 3. Computes a spatial similarity matrix `S` where
#'    \eqn{S_{ij} = exp(-\gamma * d_{ij}^2)} for the Euclidean distance
#'    between positions. Small values below `S_th` are set to zero.
#' 4. Normalizes `S` so that row sums equal 1.
#' 5. Optimizes the parameter `gamma` to minimize \eqn{(mean(diag(S)) - \tau)^2}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{counts}{Filtered and converted numeric matrix of counts.}
#'   \item{S}{Normalized spatial similarity matrix.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose counts is a numeric matrix with column names like "10x15"
#' result <- load_data(counts, filter_th = 5, tau = 0.4)
#' counts_filtered <- result$counts
#' S_matrix <- result$S
#' }
#'
#' @export
load_data <- function(counts, filter_th=10, tau=0.5, S_th=1e-3){

    filter <- rowSums(counts) > filter_th
    counts <- counts[filter, ]
    counts <- as.matrix(counts)

    positions <- matrix(as.numeric(unlist(strsplit(colnames(counts), "x"))), 
                    ncol=2, byrow = TRUE)
    x <- positions[,1]
    y <- positions[,2]

    D2 <- as.matrix(dist(cbind(x, y)))^2

    meanValue <- function(gamma, D2, tau) {
        S <- exp(-gamma * D2)
        S[S < S_th] <- 0 
        rs <- rowSums(S)
        rs[rs == 0] <- 1
        S <- Matrix::Diagonal(x = 1/rs) %*% S
        return((mean(Matrix::diag(S)) - tau)^2)
    }

    gamma <- optim(1, meanValue, method="BFGS", tau=tau, D2=D2)$par

    S <- exp(-gamma * as.matrix(dist(cbind(x,y)))^2)
    S[S < S_th] <- 0 

    S <- Matrix::Diagonal(x = 1/rowSums(S)) %*% S

    return(list(counts=counts, S=S))

}