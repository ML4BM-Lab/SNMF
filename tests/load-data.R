library(SNMF)

expect_error <- function(expr, pattern) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(tolower(pattern), tolower(conditionMessage(error)), fixed = TRUE))
}
counts <- matrix(c(1, 0, 4, 2, 0, 5, 3, 0, 6), nrow = 3,
                 dimnames = list(c("a", "zero", "b"), c("0x0", "1x0", "0x1")))
result <- load_data(counts, filter_th = 0)
stopifnot(identical(result$counts, counts[c(1, 3), , drop = FALSE]))
stopifnot(identical(dim(result$S), c(3L, 3L)))
stopifnot(all(is.finite(as.matrix(result$S))), all(as.matrix(result$S) >= 0))
stopifnot(isTRUE(all.equal(as.numeric(Matrix::rowSums(result$S)), rep(1, 3))))
# Supported storage types must preserve the output and ordering.
for (input in list(as.data.frame(counts), Matrix::Matrix(counts, sparse = TRUE))) {
    stopifnot(isTRUE(all.equal(load_data(input, filter_th = 0), result)))
}
# Filtering must retain matrix dimensions, including singleton cases.
one_feature <- load_data(counts, filter_th = 10)
stopifnot(identical(dim(one_feature$counts), c(1L, 3L)))
one_spot <- load_data(counts[, 1, drop = FALSE], filter_th = 0)
stopifnot(identical(dim(one_spot$counts), c(2L, 1L)))
stopifnot(identical(as.numeric(one_spot$S), 1))
expect_error(load_data(counts, filter_th = 100), "No features remain")

for (input in list(1:3, list(a = 1), matrix(numeric(), 0, 2), matrix(numeric(), 2, 0))) {
    expect_error(load_data(input), "nonempty numeric")
}
for (value in list(NA_real_, NaN, Inf, -1, "bad", 1i)) {
    bad <- counts
    bad[1, 1] <- value
    expect_error(load_data(bad), "finite, nonnegative numeric")
}
for (name in c("filter_th", "tau", "S_th")) {
    for (value in list(NULL, numeric(), c(0.1, 0.2), NA_real_, NaN, Inf, -1, "0.5", TRUE, 1i)) {
        args <- list(counts = counts, filter_th = 0)
        args[name] <- list(value)
        expect_error(do.call(load_data, args), paste0(name, " must be"))
    }
}
expect_error(load_data(counts, tau = 0), "tau must be")
expect_error(load_data(counts, tau = 1.1), "tau must be")
expect_error(load_data(counts, S_th = 1.1), "S_th must be")

bad <- counts
colnames(bad) <- NULL
expect_error(load_data(bad), "column names encoding coordinates")
for (name in c(NA_character_, "", "1", "x1", "1x", "1x2x", "1x2x3", "ax2", "NAx1", "Infx2")) {
    bad <- counts
    colnames(bad)[1] <- name
    expect_error(load_data(bad), if (is.na(name)) "column names" else "coordinates")
}
# Signed, decimal, and scientific notation retain the existing numeric parser.
valid <- counts
colnames(valid) <- c("-1.5x0", "1e0x0", "0x+1.5")
stopifnot(all(is.finite(as.matrix(load_data(valid, filter_th = 0)$S))))
bad <- counts
colnames(bad)[1] <- "1e200x0"
expect_error(load_data(bad, filter_th = 0), "nonfinite squared distances")
