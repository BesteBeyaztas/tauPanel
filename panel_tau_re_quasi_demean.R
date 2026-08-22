panel_tau_re_quasi_demean <- function(
    Y,
    X,
    theta,
    residual_matrix,
    tuning_constant = 4.685
) {
  N <- nrow(Y)
  TT <- ncol(Y)
  K <- dim(X)[3]
  
  if (!identical(dim(residual_matrix), dim(Y))) {
    stop("residual_matrix must have the same dimensions as Y.")
  }
  
  if (length(theta) != 1 ||
      !is.finite(theta) ||
      theta < 0 ||
      theta > 1) {
    stop("theta must be a finite scalar in [0, 1].")
  }
  
  residual_center <- apply(
    residual_matrix,
    1,
    median,
    na.rm = TRUE
  )
  
  centered_residuals <- sweep(
    residual_matrix,
    1,
    residual_center,
    "-"
  )
  
  residual_scale <- Qn(
    as.vector(centered_residuals)
  )
  
  if (!is.finite(residual_scale) || residual_scale <= 1e-8) {
    residual_scale <- mad(
      as.vector(centered_residuals),
      constant = 1.4826,
      na.rm = TRUE
    )
  }
  
  if (!is.finite(residual_scale) || residual_scale <= 1e-8) {
    residual_scale <- sqrt(
      mean(centered_residuals^2, na.rm = TRUE)
    )
  }
  
  if (!is.finite(residual_scale) || residual_scale <= 1e-8) {
    residual_scale <- 1
  }
  
  standardized_residuals <-
    centered_residuals / (tuning_constant * residual_scale)
  
  weights <- matrix(
    0,
    nrow = N,
    ncol = TT
  )
  
  inside <- abs(standardized_residuals) < 1
  weights[inside] <-
    (1 - standardized_residuals[inside]^2)^2
  
  weights[!is.finite(weights)] <- 0
  
  for (i in seq_len(N)) {
    positive <- which(weights[i, ] > 0)
    
    if (length(positive) < 2) {
      retained <- order(
        abs(centered_residuals[i, ]),
        na.last = NA
      )[seq_len(min(2, TT))]
      
      weights[i, ] <- 0
      weights[i, retained] <- 1
    }
  }
  
  Yq <- matrix(
    NA_real_,
    nrow = N,
    ncol = TT
  )
  
  Xq_arr <- array(
    NA_real_,
    dim = c(N, TT, K)
  )
  
  weighted_y_means <- numeric(N)
  weighted_x_means <- matrix(
    NA_real_,
    nrow = N,
    ncol = K
  )
  
  for (i in seq_len(N)) {
    weight_sum <- sum(weights[i, ])
    
    if (!is.finite(weight_sum) || weight_sum <= 0) {
      stop("The robust quasi-demeaning weights are degenerate.")
    }
    
    weighted_y_means[i] <-
      sum(weights[i, ] * Y[i, ]) / weight_sum
    
    Yq[i, ] <-
      Y[i, ] - theta * weighted_y_means[i]
    
    for (k in seq_len(K)) {
      weighted_x_means[i, k] <-
        sum(weights[i, ] * X[i, , k]) / weight_sum
      
      Xq_arr[i, , k] <-
        X[i, , k] - theta * weighted_x_means[i, k]
    }
  }
  
  Xq <- matrix(
    NA_real_,
    nrow = N * TT,
    ncol = K
  )
  
  for (k in seq_len(K)) {
    Xq[, k] <- c(t(Xq_arr[, , k]))
  }
  
  colnames(Xq) <- paste0("X", seq_len(K), "_tau_re_quasi")
  
  list(
    y = c(t(Yq)),
    x = Xq,
    id = rep(seq_len(N), each = TT),
    theta = theta,
    weights = weights,
    residual_scale = residual_scale,
    weighted_y_means = weighted_y_means,
    weighted_x_means = weighted_x_means
  )
}