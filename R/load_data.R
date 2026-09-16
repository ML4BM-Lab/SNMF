
#' Load and preprocess count data for NMF
#'
#' This function filters and preprocesses a count matrix for non-negative
#' matrix factorization (NMF), and computes a spatial similarity matrix `S`
#' based on the positions of the columns.
#'
#' @param counts Nonempty numeric matrix, data frame, or numeric Matrix object
#'   containing finite, nonnegative counts, where rows are features (e.g., genes)
#'   and columns are observations (e.g., spots or cells). Column names should
#'   encode positions in the form "<coord>x<coord>" (for example, "10x15").
#' @param filter_th Finite nonnegative number, minimum sum per row to keep a feature (default 10).
#'   Features with row sums ≤ `filter_th` are removed.
#' @param tau Finite number in (0, 1], target mean value for optimization of gamma (default 0.8).
#'   This is used in computing the spatial similarity matrix.
#' @param S_th Finite number in [0, 1], threshold below which entries of the similarity matrix
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
load_data <- function(counts, filter_th=10, tau=0.8, S_th=1e-3){

    if (!(is.matrix(counts) || is.data.frame(counts) || inherits(counts, "Matrix")) ||
        length(dim(counts)) != 2L || any(dim(counts) == 0L)) {
        stop("counts must be a nonempty numeric matrix, data frame, or Matrix object.", call. = FALSE)
    }
    values <- as.matrix(counts)
    if (!is.numeric(values) || is.complex(values) ||
        any(!is.finite(values)) || any(values < 0)) {
        stop("counts must contain only finite, nonnegative numeric values.", call. = FALSE)
    }

    rm(values)

    check_scalar <- function(value, name, lower, upper = Inf, open_lower = FALSE) {
        if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
            !is.finite(value) || value < lower || value > upper ||
            (open_lower && value == lower)) {
            interval <- if (open_lower) "(" else "["
            stop(name, " must be a finite numeric scalar in ", interval,
                 lower, ", ", upper, "].", call. = FALSE)
        }
    }
    check_scalar(filter_th, "filter_th", 0)
    check_scalar(tau, "tau", 0, 1, open_lower = TRUE)
    check_scalar(S_th, "S_th", 0, 1)

    spot_names <- colnames(counts)
    if (is.null(spot_names) || anyNA(spot_names)) {
        stop('counts must have column names encoding coordinates, for example "10x15".', call. = FALSE)
    }
    coordinate_parts <- strsplit(spot_names, "x", fixed = TRUE)
    # Reject trailing separators too: strsplit() drops trailing empty fields.
    valid_parts <- vapply(coordinate_parts, length, integer(1)) == 2L
    if (any(!valid_parts) || any(endsWith(spot_names, "x"))) {
        stop('Each counts column name must contain exactly two numeric coordinates separated by "x".', call. = FALSE)
    }
    coordinates <- suppressWarnings(as.numeric(unlist(coordinate_parts)))
    if (any(!is.finite(coordinates))) {
        stop("Coordinates in counts column names must be finite numbers.", call. = FALSE)
    }

    filter <- rowSums(counts) > filter_th
    if (!any(filter)) {
        stop("No features remain after filtering; reduce filter_th or check counts.", call. = FALSE)
    }
    counts <- counts[filter, , drop = FALSE]
    counts <- as.matrix(counts)

    positions <- matrix(coordinates, ncol=2, byrow = TRUE)
    x <- positions[,1]
    y <- positions[,2]

    D2 <- as.matrix(stats::dist(cbind(x, y)))^2
    if (any(!is.finite(D2))) {
        stop("Coordinates produce nonfinite squared distances; check their magnitude.", call. = FALSE)
    }

    meanValue <- function(gamma, D2, tau) {
        S <- exp(-gamma * D2)
        S[S < S_th] <- 0 
        rs <- rowSums(S)
        rs[rs == 0] <- 1
        S <- Matrix::Diagonal(x = 1/rs) %*% S
        return((mean(Matrix::diag(S)) - tau)^2)
    }

    gamma <- stats::optim(1, meanValue, method="L-BFGS-B", lower=1e-12, tau=tau, D2=D2)$par

    S <- exp(-gamma * as.matrix(stats::dist(cbind(x,y)))^2)
    S[S < S_th] <- 0 

    if (any(!is.finite(S)) || any(rowSums(S) <= 0)) {
        stop("Spatial similarity calculation produced nonfinite values or empty rows.", call. = FALSE)
    }
    S <- Matrix::Diagonal(x = 1/rowSums(S)) %*% S

    return(list(counts=counts, S=S))

}