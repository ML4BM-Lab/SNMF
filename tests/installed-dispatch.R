# Run from a fresh R session against the installed package.
# GPU integration is opt-in because it requires a configured torch/CUDA backend:
# SNMF_TEST_GPU=true Rscript --vanilla tests/installed-dispatch.R
library(SNMF)
stopifnot(!"package:GPUmatrix" %in% search())

ns <- asNamespace("SNMF")
for (name in c("as.matrix", "colSums", "mean", "ncol", "nrow", "rowSums", "t")) {
    stopifnot(identical(get(name, envir = ns), getExportedValue("GPUmatrix", name)))
}

# CPU reference for the NB multiplicative updates, independent of GPU dispatch.
V <- matrix(seq_len(12), 3, 4)
W <- matrix(seq_len(6) / 6, 3, 2)
H <- matrix(seq_len(8) / 8, 2, 4)
S <- diag(4)
WHS <- W %*% H %*% S
R <- list(phi = matrix(0.1, 3, 4), eps = 1e-10, phi.min = 1e-8)
ratio <- V / WHS
q <- (1 + R$phi * V) / (1 + R$phi * WHS)
expected_W <- W * (ratio %*% t(H %*% S)) / (q %*% t(H %*% S))
expected_H <- H * (t(W) %*% ratio %*% t(S)) / (t(W) %*% q %*% t(S))
stopifnot(isTRUE(all.equal(SNMF:::updateW(V, W, H, S, WHS, R), expected_W)))
stopifnot(isTRUE(all.equal(SNMF:::updateH(V, W, H, S, WHS, R), expected_H)))

if (identical(Sys.getenv("SNMF_TEST_GPU"), "true")) {
    data("tnbc", package = "SNMF", envir = environment())
    input <- load_data(tnbc[, seq_len(20), drop = FALSE])
    run <- function(fun) {
        # Exercise random-start selection, updateR(), and iteration-10 rescaling.
        suppressWarnings(fun(input$counts, input$S, k = 2, niter = 20,
                             Rupdate_iter = 2, num_initializations = 2,
                             tol = 0, seed = 42))
    }
    installed <- run(SNMF::snmf)
    rng_installed <- .Random.seed
    stopifnot(!"package:GPUmatrix" %in% search())
    stopifnot(all(vapply(installed, function(x) all(is.finite(x)), logical(1))))
    stopifnot(identical(dim(installed$W), c(nrow(input$counts), 2L)))
    stopifnot(identical(dim(installed$H), c(2L, ncol(input$counts))))

    # Recreate the working script's lookup environment with GPUmatrix attached.
    # Using the same function bodies isolates namespace/dispatch differences.
    library(GPUmatrix)
    script <- new.env(parent = globalenv())
    for (name in c("snmf", "factorize", "updateW", "updateH", "updateR",
                   "safe_pmax", "controlDimensionNMF")) {
        fun <- get(name, envir = ns)
        environment(fun) <- script
        assign(name, fun, envir = script)
    }
    reference <- run(script$snmf)
    stopifnot(identical(.Random.seed, rng_installed))
    stopifnot(isTRUE(all.equal(installed, reference, tolerance = 1e-6)))
}
