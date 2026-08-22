apply_panel_contamination <- function(X, epsilon,
                                      scheme = c("none",
                                                 "random_vertical",
                                                 "block_vertical",
                                                 "random_leverage",
                                                 "block_leverage",
                                                 "random_vertical_leverage"),
                                      cont = 0,
                                      eps_mean = 10,
                                      eps_sd = 1,
                                      x_mean = 5,
                                      x_sd = 1,
                                      block_fraction = 0.5) {
  
  scheme <- match.arg(scheme)
  
  N_ <- nrow(epsilon)
  T_ <- ncol(epsilon)
  K  <- dim(X)[3]
  
  if (scheme == "none" || cont <= 0) {
    
    return(list(X = X,
                epsilon = epsilon,
                outlier_cells = integer(0),
                outlier_units = integer(0)))
  }
  
  m_target <- max(1L, round(cont * N_ * T_))
  
  if (scheme %in%
      c("random_vertical", "random_leverage",
        "random_vertical_leverage")) {
    
    cells <- sample.int(N_ * T_,
                        size = m_target,
                        replace = FALSE)
    
    unit_index <- ((cells - 1) %/% T_) + 1
    time_index <- ((cells - 1) %% T_) + 1
    
  } else {
    
    h <- max(1,
             ceiling(block_fraction * T_))
    
    n_units <- min(N_,
                   max(1, round(m_target / h)))
    
    selected_units <- sample.int(N_,
                                 size = n_units,
                                 replace = FALSE)
    
    block_pairs <- lapply(selected_units,
                          function(i) {
                            
                            # Consecutive block of h time points.
                            start_time <- sample.int(T_ - h + 1,
                                                     size = 1)
                            
                            times <- start_time:(start_time + h - 1)
                            
                            cbind(unit = rep(i, h),
                                  time = times)
                          }
    )
    
    block_pairs <- do.call(rbind, block_pairs)
    
    unit_index <- block_pairs[, "unit"]
    time_index <- block_pairs[, "time"]
    
    cells <- (unit_index - 1L) * T_ + time_index
  }
  
  
  if (scheme %in%
      c("random_vertical", "block_vertical",
        "random_vertical_leverage")) {
    
    epsilon_vec <- c(t(epsilon))
    
    epsilon_vec[cells] <- rnorm(length(cells),
                                mean = eps_mean,
                                sd = eps_sd)
    
    epsilon <- matrix(epsilon_vec,
                      nrow = N_,
                      ncol = T_,
                      byrow = TRUE)
  }
  
  if (scheme %in%
      c("random_leverage", "block_leverage",
        "random_vertical_leverage")) {
    
    for (j in seq_along(cells)) {
      
      X[unit_index[j],
        time_index[j],
        seq_len(K)] <- rnorm(K,
                             mean = x_mean,
                             sd = x_sd)
    }
  }
  
  return(list(X = X,
              epsilon = epsilon,
              outlier_cells = sort(unique(cells)),
              outlier_units = sort(unique(unit_index))))
}