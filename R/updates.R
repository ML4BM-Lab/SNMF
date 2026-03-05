updateW <- function (V, W, H, S, R, S2) {
  # R is the reconstructed matrix: W %*% H %*% S
  # S2 is a matrix of ones with the same dimensions as V
  
  tHS <- t(H %*% S) # Pre-calculate t(H %*% S) for efficiency
  num <- (V / R) %*% tHS # Numerator of the update rule
  denom <- S2 %*% tHS # Denominator of the update rule
  W <- W * num / denom # Multiplicative update rule for W
  return(W)
}

# Optimized version of updateH function
updateH <- function (V, W, H, S, R, S1) {
  # S1 is (a matrix of ones with dimensions nrow(W) x nrow(S)) %*% S
  
  num <- (t(W) %*% (V / R)) %*% t(S) # Numerator of the update rule
  denom <- t(W) %*% S1 # Denominator of the update rule
  H <- H * num / denom # Multiplicative update rule for H
  return(H)
}