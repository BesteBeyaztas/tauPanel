panel_within_transform <- function(Y, X, robust_center = FALSE) {
  N <- nrow(Y)
  TT <- ncol(Y)
  K <- dim(X)[3]
  
  Yw <- matrix(NA, nrow = N, ncol = TT)
  Xw_arr <- array(NA, dim = c(N, TT, K))
  
  center_fun <- function(z) {
    if (robust_center) median(z, na.rm = TRUE) else mean(z, na.rm = TRUE)
  }
  
  for (i in 1:N) {
    Yw[i, ] <- Y[i, ] - center_fun(Y[i, ])
    
    for (k in 1:K) {
      Xw_arr[i, , k] <- X[i, , k] - center_fun(X[i, , k])
    }
  }
  
  Xw <- matrix(NA, nrow = N * TT, ncol = K)
  for (k in 1:K) {
    Xw[, k] <- c(t(Xw_arr[, , k]))
  }
  
  colnames(Xw) <- paste0("X", 1:K, "_within")
  
  list(
    y = c(t(Yw)),
    x = Xw,
    id = rep(seq_len(N), each = TT)
  )
}