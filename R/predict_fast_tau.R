predict_fast_tau <- function(object, newdata) {
  
  if (!is.list(object) || is.null(object$beta)) {
    stop(
      "object must be a fitted fast tau object containing an estimated ",
      "coefficient vector named 'beta'."
    )
  }
  
  predict_tau_structural_mean(
    beta_hat = object$beta,
    newdata = newdata
  )
}
