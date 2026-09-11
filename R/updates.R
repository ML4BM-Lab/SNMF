updateW <- function(V, W, H, S, WHS, R) {

  eps <- if (!is.null(R$eps)) R$eps else 1e-10

  WHS <- safe_pmax(WHS, eps)

  Phi <- R$phi
  Phi <- safe_pmax(Phi, if (!is.null(R$phi.min)) R$phi.min else 1e-8)

  Rpos <- V / WHS

  Q <- (1 + Phi * V) / (1 + Phi * WHS)

  B <- H %*% S

  numerator <- Rpos %*% t(B)
  denominator <- Q %*% t(B)

  W.new <- W * numerator / safe_pmax(denominator, eps)
  W.new <- safe_pmax(W.new, eps)

  return(W.new)
}

# Optimized version of updateH function
updateH <- function(V, W, H, S, WHS, R) {

  eps <- if (!is.null(R$eps)) R$eps else 1e-10

  WHS <- safe_pmax(WHS, eps)

  Phi <- R$phi
  Phi <- safe_pmax(Phi, if (!is.null(R$phi.min)) R$phi.min else 1e-8)

  Rpos <- V / WHS

  Q <- (1 + Phi * V) / (1 + Phi * WHS)

  numerator   <- t(W) %*% Rpos %*% t(S)
  denominator <- t(W) %*% Q %*% t(S)

  H.new <- H * numerator / safe_pmax(denominator, eps)
  H.new <- safe_pmax(H.new, eps)

  return(H.new)
}

updateR <- function(
  V, W, H, S, WHS, R
) {

  eps <- if (!is.null(R$eps)) R$eps else 1e-10
  phi.min <- if (!is.null(R$phi.min)) R$phi.min else 1e-8
  phi.max <- if (!is.null(R$phi.max)) R$phi.max else 1e3
  inner.iter <- if (!is.null(R$inner.iter)) R$inner.iter else 10

  m <- nrow(V)  # genes
  n <- ncol(V)  # spots

  WHS <- safe_pmax(WHS, eps)

  # Moment target:
  # Var(V_ij) = mu_ij + phi_ij * mu_ij^2
  D <- (V - WHS)^2 - WHS
  C <- WHS^2

  if (is.null(R$alpha)) {
    R$alpha <- rep(0, m)
  }

  if (is.null(R$beta)) {
    R$beta <- rep(0, n)
  }

  alpha <- R$alpha
  beta  <- R$beta

  for (iter in seq_len(inner.iter)) {

    alpha <- (
      rowSums(D) -
      rowSums(C * matrix(beta, nrow = m, ncol = n, byrow = TRUE))
    ) / safe_pmax(rowSums(C), eps)

    beta <- (
      colSums(D) -
      colSums(C * matrix(alpha, nrow = m, ncol = n, byrow = FALSE))
    ) / safe_pmax(colSums(C), eps)

    # Identifiability constraint: mean(beta) = 0
    beta.mean <- mean(beta)

    beta  <- beta - beta.mean
    alpha <- alpha + beta.mean
  }

  Phi <- outer(alpha, beta, "+")

  Phi <- pmin(safe_pmax(Phi, phi.min), phi.max)

  R$alpha <- alpha
  R$beta  <- beta
  R$phi   <- Phi
  R$theta <- 1 / Phi

  return(R)
  
}