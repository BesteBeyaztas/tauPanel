fast_tau_re <- function(
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
  validate_panel_data(Y, X)
  
  result <- fit_tau_re_estimator(
    Y = Y,
    X = X,
    N = N,
    kk = kk,
    tt = tt,
    rr = rr,
    approximate = approximate,
    cluster_start = cluster_start,
    m_cluster = m_cluster,
    seed = seed
  )
  
  result$effect <- "random"
  class(result) <- c("fast_tau_re_fit", "fast_tau_panel_fit")
  result
}