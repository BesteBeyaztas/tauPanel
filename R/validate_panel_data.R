validate_panel_data <- function(Y, X) {
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric N by T matrix.")
  }
  
  if (length(dim(X)) != 3 || !is.numeric(X)) {
    stop("X must be a numeric N by T by K array.")
  }
  
  if (dim(X)[1] != nrow(Y) || dim(X)[2] != ncol(Y)) {
    stop("The first two dimensions of X must match the dimensions of Y.")
  }
  
  if (nrow(Y) < 2 || ncol(Y) < 2 || dim(X)[3] < 1) {
    stop("The panel must contain at least two individuals, two time points, and one regressor.")
  }
  
  if (any(!is.finite(Y)) || any(!is.finite(X))) {
    stop("Y and X must contain only finite values.")
  }
  
  invisible(TRUE)
}
