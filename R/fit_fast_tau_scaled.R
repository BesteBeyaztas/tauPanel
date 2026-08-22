fit_fast_tau_scaled <- function(x, y, ...) {
  x <- as.matrix(x)
  y <- as.vector(y)
  
  if (nrow(x) != length(y)) {
    stop("The number of rows in x must equal the length of y.")
  }
  
  x_scale <- apply(
    x,
    2,
    function(z) {
      mad(
        z,
        constant = 1.4826,
        na.rm = TRUE
      )
    }
  )
  
  fallback_scale <- sqrt(
    colMeans(x^2, na.rm = TRUE)
  )
  
  invalid_scale <- !is.finite(x_scale) |
    x_scale < 1e-8
  
  x_scale[invalid_scale] <-
    fallback_scale[invalid_scale]
  
  x_scale[
    !is.finite(x_scale) |
      x_scale < 1e-8
  ] <- 1
  
  x_scaled <- sweep(
    x,
    2,
    x_scale,
    "/"
  )
  
  fit <- FastTau(
    x = x_scaled,
    y = y,
    ...
  )
  
  beta_scaled <- as.numeric(fit$beta)
  
  fit$beta_scaled <- beta_scaled
  fit$x_scale <- x_scale
  fit$beta <- beta_scaled / x_scale
  
  fit
}
