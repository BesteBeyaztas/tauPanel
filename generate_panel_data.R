generate_panel_data <- function(
    N_ = 100,
    T_ = 3,
    K = 4,
    distribution = "normal",
    beta = c(-3.4, 1.5, 2.8, -1.2),
    sigma2_a = 1,
    sigma2_e = 1,
    effect_model = c("re", "fe"),
    alpha_corr = 0.6,
    gamma_alpha = NULL,
    
    contam_scheme = c(
      "none",
      "random_vertical",
      "block_vertical",
      "random_leverage",
      "block_leverage",
      "random_vertical_leverage"
    ),
    
    eps_mean = 10,
    eps_std = 1,
    cont = 0.05,
    mu_leverage = 5,
    std_leverage = 1,
    block_fraction = 0.5,
    chi_df = NULL
) {
  
  effect_model_all <- c("re", "fe")
  effect_model <- effect_model_all[pmatch(tolower(effect_model[1]),
                                          effect_model_all)]
  
  if (length(effect_model) != 1 || is.na(effect_model)) {
    stop("effect_model must be either 're' or 'fe'.")
  }
  
  if (length(beta) != K || any(!is.finite(beta))) {
    stop("beta must be a finite numeric vector of length K.")
  }
  
  if (is.null(chi_df)) {
    chi_df <- if (effect_model == "fe") 2 else 1
  }
  
  if (length(chi_df) != 1 || !is.finite(chi_df) || chi_df <= 0) {
    stop("chi_df must be a positive finite scalar.")
  }
  
  contam_scheme <- match.arg(
    tolower(contam_scheme[1]),
    choices = c(
      "none",
      "random_vertical",
      "block_vertical",
      "random_leverage",
      "block_leverage",
      "random_vertical_leverage"
    )
  )
  
  X = array(rnorm(N_ * T_ * K), dim = c(N_, T_, K))
  
  X[, , 1] = matrix(rchisq(N_ * T_, df = 2) - 2, nrow = N_)
  
  X_clean <- X
  
  if (effect_model == "re") {
    
    alpha = matrix(rnorm(N_, sd = sqrt(sigma2_a)), nrow = N_)
    
  } else if (effect_model == "fe") {
    
    X_bar <- matrix(
      NA_real_,
      nrow = N_,
      ncol = K
    )
    
    for (k in seq_len(K)) {
      X_bar[, k] <- rowMeans(X[, , k])
    }
    
    if (is.null(gamma_alpha)) {
      
      if (K == 4) {
        gamma_alpha <- c(1.0, -0.7, 0.5, 0.3)
      } else {
        gamma_alpha <- seq(1, 0.4, length.out = K)
      }
      
    } else {
      
      if (length(gamma_alpha) != K) {
        stop("gamma_alpha must have length K.")
      }
      
      if (any(!is.finite(gamma_alpha))) {
        stop("gamma_alpha must contain only finite values.")
      }
    }
    
    gamma_norm <- sqrt(sum(gamma_alpha^2))
    
    if (!is.finite(gamma_norm) || gamma_norm <= 0) {
      stop("gamma_alpha must not be the zero vector.")
    }
    
    gamma_alpha <- gamma_alpha / gamma_norm
    
    z <- as.vector(X_bar %*% gamma_alpha)
    z <- as.numeric(scale(z))
    
    if (any(!is.finite(z))) {
      stop("The Mundlak index z is degenerate. Check N_, T_, K, and gamma_alpha.")
    }
    
    eta <- rnorm(N_)
    
    eta <- eta -
      as.numeric(
        crossprod(z, eta) /
          crossprod(z, z)
      ) * z
    
    eta <- as.numeric(scale(eta))
    
    if (any(!is.finite(eta))) {
      stop("The orthogonalized eta vector is degenerate. Increase N_ or check the FE DGP.")
    }
    
    rho <- max(
      min(alpha_corr, 0.99),
      -0.99
    )
    
    alpha_vec <- sqrt(sigma2_a) * (
      rho * z +
        sqrt(1 - rho^2) * eta
    )
    
    alpha <- matrix(
      alpha_vec,
      nrow = N_
    )
  }
  
  Y  <- matrix(NA_real_, nrow = N_, ncol = T_)
  MU <- matrix(NA_real_, nrow = N_, ncol = T_)
  
  epsilon <- generate_epsilon(
    N_ = N_,
    T_ = T_,
    distribution = distribution,
    sigma2_e = sigma2_e,
    chi_df = chi_df
  )
  
  contam <- apply_panel_contamination(
    X = X,
    epsilon = epsilon,
    scheme = contam_scheme,
    cont = cont,
    eps_mean = eps_mean,
    eps_sd = eps_std,
    x_mean = mu_leverage,
    x_sd = std_leverage,
    block_fraction = block_fraction
  )
  
  X <- contam$X
  epsilon <- contam$epsilon
  
  outlier_cells <- contam$outlier_cells
  outlier_units <- contam$outlier_units
  
  if (contam_scheme %in% c(
    "random_leverage",
    "block_leverage",
    "random_vertical_leverage"
  )) {
    X_for_y <- X_clean
  } else {
    X_for_y <- X
  }
  
  for (t in seq_len(T_)) {
    MU[, t] <- as.vector(X_for_y[, t, ] %*% beta)
    Y[, t]  <- MU[, t] + as.vector(alpha) + epsilon[, t]
  }
  
  return(list(
    X = X,
    Y = Y,
    MU = MU,
    alpha = alpha,
    epsilon = epsilon,
    outlier_cells = outlier_cells,
    outlier_units = outlier_units,
    gamma_alpha = gamma_alpha,
    effect_model = effect_model,
    chi_df = chi_df
  ))
}