fast_tau_fe <- function(
    Y,
    X,
    transform = c("first_difference", "within_mean", "within_median"),
    N = 2000,
    kk = 2,
    tt = 20,
    rr = 2,
    approximate = 0,
    cluster_start = TRUE,
    m_cluster = NULL,
    seed = NULL
) {
  validate_panel_data(Y, X)
  transform <- match.arg(transform)
  
  transformed <- switch(
    transform,
    first_difference = panel_first_difference(Y, X),
    within_mean = panel_within_transform(Y, X, robust_center = FALSE),
    within_median = panel_within_transform(Y, X, robust_center = TRUE)
  )
  
  result <- fit_tau_fe_estimator(
    y = transformed$y,
    x = transformed$x,
    id = transformed$id,
    N = N,
    kk = kk,
    tt = tt,
    rr = rr,
    approximate = approximate,
    cluster_start = cluster_start,
    m_cluster = m_cluster,
    seed = seed
  )
  
  result$effect <- "fixed"
  result$transform <- transform
  result$transformed_data <- transformed
  class(result) <- c("fast_tau_fe_fit", "fast_tau_panel_fit")
  result
}
