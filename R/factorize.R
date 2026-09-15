factorize <- function(
  V,
  S = diag(ncol(V)),
  k = 10,
  Winit = NULL,
  Hinit = NULL,
  tol = 1e-03,
  Rupdate_iter = 20,
  niter = 100,
  num_initializations = 10
)
{

  dtype <- "float32"

  objectClass <- class(V)[[1]]
  objectPackage <- attr(class(V), "package")

  if (!is.null(objectPackage) &&
      (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {

    if (is.null(Winit)) {
      Winit <- GPUmatrix::gpu.matrix(
        stats::runif(nrow(V) * k), nrow(V), k,
        dtype = GPUmatrix::dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )
    }

    if (is.null(Hinit)) {
      Hinit <- GPUmatrix::gpu.matrix(
        stats::runif(k * ncol(V)), k, ncol(V),
        dtype = GPUmatrix::dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )
    }

  } else {

    if (is.null(Winit)) {
      Winit <- matrix(stats::runif(nrow(V) * k), nrow(V), k)
    }

    if (is.null(Hinit)) {
      Hinit <- matrix(stats::runif(k * ncol(V)), k, ncol(V))
    }
  }

  controlDimensionNMF(Winit, Hinit, V, k)

  Vold <- V

  if (!is.null(objectPackage) &&
      (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {

    R <- list(
      alpha = rep(0, nrow(V)),
      beta = rep(0, ncol(V)),
      phi = GPUmatrix::gpu.matrix(
        matrix(1e-4, nrow(V), ncol(V)),
        dtype = GPUmatrix::dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      ),
      eps = 1e-10,
      phi.min = 1e-8,
      phi.max = 1e3,
      inner.iter = 10
    )

  } else {

    R <- list(
      alpha = rep(0, nrow(V)),
      beta = rep(0, ncol(V)),
      phi = matrix(1e-4, nrow(V), ncol(V)),
      eps = 1e-10,
      phi.min = 1e-8,
      phi.max = 1e3,
      inner.iter = 10
    )
  }

  initial_iterations <- max(1, floor(niter / 10))

  best_loss <- Inf
  best_W <- NULL
  best_H <- NULL

  for (init_run in seq_len(num_initializations)) {

    if (!is.null(objectPackage) &&
        (objectClass == "gpu.matrix.torch" || objectClass == "gpu.matrix.tensorflow")) {

      W_current <- GPUmatrix::gpu.matrix(
        stats::runif(nrow(V) * k), nrow(V), k,
        dtype = GPUmatrix::dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )

      H_current <- GPUmatrix::gpu.matrix(
        stats::runif(k * ncol(V)), k, ncol(V),
        dtype = GPUmatrix::dtype(V),
        type = GPUmatrix:::typeGPUmatrix(V),
        device = GPUmatrix:::device(V)
      )

    } else {

      W_current <- matrix(stats::runif(nrow(V) * k), nrow(V), k)
      H_current <- matrix(stats::runif(k * ncol(V)), k, ncol(V))
    }

    for (iter_init in seq_len(initial_iterations)) {

      WHS_current <- W_current %*% H_current %*% S
      WHS_current <- safe_pmax(WHS_current, R$eps)

      W_current <- updateW(V, W_current, H_current, S, WHS_current, R)

      WHS_current <- W_current %*% H_current %*% S
      WHS_current <- safe_pmax(WHS_current, R$eps)

      H_current <- updateH(V, W_current, H_current, S, WHS_current, R)
    }

    V_reconstructed <- W_current %*% H_current %*% S
    current_loss <- mean((V_reconstructed - V)^2)

    if (current_loss < best_loss) {
      best_loss <- current_loss
      best_W <- W_current
      best_H <- H_current
    }
  }

  Winit <- best_W
  Hinit <- best_H

  for (iter in seq_len(niter)) {

    WHS <- Winit %*% Hinit %*% S
    WHS <- safe_pmax(WHS, R$eps)

    Winit <- updateW(V, Winit, Hinit, S, WHS, R)

    WHS <- Winit %*% Hinit %*% S
    WHS <- safe_pmax(WHS, R$eps)

    Hinit <- updateH(V, Winit, Hinit, S, WHS, R)

    if (iter %% Rupdate_iter == 0) {
      WHS <- Winit %*% Hinit %*% S
      WHS <- safe_pmax(WHS, R$eps)

      R <- updateR(V, Winit, Hinit, S, WHS, R)
    }

    if (iter %% 100 == 0) {
      cat("Iteration:", iter, "\n")
    }

    if (iter %% 10 == 0) {

      myD <- colSums(Winit)
      myD <- safe_pmax(myD, R$eps)

      Winit <- t(t(Winit) / myD)
      Hinit <- Hinit * myD

      Vnew <- Winit %*% Hinit %*% S

      mse_change <- mean((Vnew - Vold)^2)

      if (is.na(mse_change)) {
        stop("Error in calculating mean squared error. Check dimensions of Vnew and Vold.")
      }

      if (mse_change < tol) {
        message("NMF converged early.")
        return(list(W = Winit, H = Hinit, phi = R$phi, alpha = R$alpha, beta = R$beta, niter=iter + num_initializations*initial_iterations))
      }

      Vold <- Vnew
    }
  }

  warning("Maximum number of iterations reached without convergence. Consider increasing 'niter'.")
  
  return(list(W = Winit, H = Hinit, phi = R$phi, alpha = R$alpha, beta = R$beta, niter = iter + num_initializations*initial_iterations))
}