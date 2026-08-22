predict_tau_structural_mean <- function(beta_hat, newdata) {
  
  beta_hat <- as.numeric(beta_hat)
  K <- length(beta_hat)
  
  if (K < 1 || any(!is.finite(beta_hat))) {
    stop("The fitted coefficient vector is invalid.")
  }
  
  if (is.atomic(newdata) && is.null(dim(newdata))) {
    
    if (!is.numeric(newdata)) {
      stop("newdata must be numeric.")
    }
    
    if (length(newdata) != K) {
      stop(
        "The predictor vector must have length ",
        K,
        ", equal to the number of fitted coefficients."
      )
    }
    
    return(as.numeric(crossprod(beta_hat, newdata)))
  }
  
  if (is.matrix(newdata) || is.data.frame(newdata)) {
    
    if (is.data.frame(newdata)) {
      numeric_columns <- vapply(newdata, is.numeric, logical(1))
      if (!all(numeric_columns)) {
        stop("All columns of newdata must be numeric.")
      }
    }
    
    X_new <- as.matrix(newdata)
    
    if (!is.numeric(X_new)) {
      stop("newdata must be numeric.")
    }
    
    if (ncol(X_new) != K) {
      stop(
        "newdata must have ",
        K,
        " columns, equal to the number of fitted coefficients."
      )
    }
    
    return(as.vector(X_new %*% beta_hat))
  }
  
  if (is.array(newdata) && length(dim(newdata)) == 3) {
    
    if (!is.numeric(newdata)) {
      stop("newdata must be a numeric array.")
    }
    
    dimensions <- dim(newdata)
    N_new <- dimensions[1]
    T_new <- dimensions[2]
    K_new <- dimensions[3]
    
    if (K_new != K) {
      stop(
        "The third dimension of newdata must be ",
        K,
        ", equal to the number of fitted coefficients."
      )
    }
    
    predicted <- matrix(
      NA_real_,
      nrow = N_new,
      ncol = T_new
    )
    
    for (t in seq_len(T_new)) {
      X_t <- matrix(
        newdata[, t, ],
        nrow = N_new,
        ncol = K
      )
      
      predicted[, t] <- as.vector(X_t %*% beta_hat)
    }
    
    new_dimnames <- dimnames(newdata)
    if (!is.null(new_dimnames)) {
      dimnames(predicted) <- new_dimnames[1:2]
    }
    
    return(predicted)
  }
  
  stop(
    "newdata must be a numeric vector, an n by K matrix/data frame, ",
    "or an N by T by K array."
  )
}
