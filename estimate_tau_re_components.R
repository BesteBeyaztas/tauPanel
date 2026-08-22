estimate_tau_re_components <- function(
    Y,
    X,
    N = 2000,
    kk = 2,
    tt = 20,
    rr = 2,
    approximate = 0,
    cluster_start = TRUE,
    m_cluster = NULL,
    seed = NULL
) {
  N_individual <- nrow(Y)
  TT <- ncol(Y)
  K <- dim(X)[3]
  
  difference_data <- panel_first_difference(
    Y = Y,
    X = X
  )
  
  preliminary_fit <- fit_fast_tau_scaled(
    x = difference_data$x,
    y = difference_data$y,
    N = N,
    kk = kk,
    tt = tt,
    rr = rr,
    approximate = approximate,
    seed = seed,
    id = difference_data$id,
    cluster_start = cluster_start,
    m_cluster = m_cluster
  )
  
  beta_initial <- as.numeric(preliminary_fit$beta)
  
  residual_matrix <- matrix(
    NA_real_,
    nrow = N_individual,
    ncol = TT
  )
  
  for (t in seq_len(TT)) {
    X_t <- matrix(
      X[, t, ],
      nrow = N_individual,
      ncol = K
    )
    
    residual_matrix[, t] <-
      Y[, t] - as.vector(X_t %*% beta_initial)
  }
  
  pair_index <- combn(seq_len(TT), 2)
  number_pairs <- ncol(pair_index)
  
  residual_differences <- numeric(N_individual * number_pairs)
  residual_sums <- numeric(N_individual * number_pairs)
  
  for (j in seq_len(number_pairs)) {
    index <- ((j - 1) * N_individual + 1):(j * N_individual)
    time_1 <- pair_index[1, j]
    time_2 <- pair_index[2, j]
    
    residual_differences[index] <-
      residual_matrix[, time_1] - residual_matrix[, time_2]
    
    residual_sums[index] <-
      residual_matrix[, time_1] + residual_matrix[, time_2]
  }
  
  residual_differences <-
    residual_differences[is.finite(residual_differences)]
  
  residual_sums <-
    residual_sums[is.finite(residual_sums)]
  
  if (length(residual_differences) < 3 ||
      length(residual_sums) < 3) {
    stop("Not enough finite residual pairs to estimate RE components.")
  }
  
  difference_scale <- Qn(residual_differences)
  sum_scale <- Qn(residual_sums)
  
  if (!is.finite(difference_scale) || difference_scale <= 0) {
    difference_scale <- mad(
      residual_differences,
      constant = 1.4826,
      na.rm = TRUE
    )
  }
  
  if (!is.finite(sum_scale) || sum_scale <= 0) {
    sum_scale <- mad(
      residual_sums,
      constant = 1.4826,
      na.rm = TRUE
    )
  }
  
  if (!is.finite(difference_scale) || difference_scale <= 0 ||
      !is.finite(sum_scale) || sum_scale <= 0) {
    stop("Robust residual scales could not be estimated.")
  }
  
  sigma2_e <- max(
    difference_scale^2 / 2,
    1e-8
  )
  
  sigma2_a <- max(
    (sum_scale^2 - difference_scale^2) / 4,
    1e-8
  )
  
  theta <- 1 - sqrt(
    sigma2_e / (sigma2_e + TT * sigma2_a)
  )
  
  if (!is.finite(theta)) {
    stop("The robust tau quasi-demeaning parameter is not finite.")
  }
  
  theta <- min(max(theta, 0), 1)
  
  list(
    beta_initial = beta_initial,
    sigma2_a = sigma2_a,
    sigma2_e = sigma2_e,
    theta = theta,
    residual_matrix = residual_matrix,
    preliminary_fit = preliminary_fit
  )
}