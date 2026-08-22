panel_first_difference <- function(Y, X) {
  N <- nrow(Y)
  TT <- ncol(Y)
  K <- dim(X)[3]
  
  if (TT < 2) {
    stop("At least two time points are required for first differences.")
  }
  
  Yd <- Y[, 2:TT, drop = FALSE] -
    Y[, 1:(TT - 1), drop = FALSE]
  
  Xd_arr <- array(
    NA_real_,
    dim = c(N, TT - 1, K)
  )
  
  for (k in seq_len(K)) {
    Xd_arr[, , k] <- X[, 2:TT, k] -
      X[, 1:(TT - 1), k]
  }
  
  Xd <- matrix(
    NA_real_,
    nrow = N * (TT - 1),
    ncol = K
  )
  
  for (k in seq_len(K)) {
    Xd[, k] <- c(t(Xd_arr[, , k]))
  }
  
  colnames(Xd) <- paste0("X", seq_len(K), "_difference")
  
  list(
    y = c(t(Yd)),
    x = Xd,
    id = rep(seq_len(N), each = TT - 1)
  )
}