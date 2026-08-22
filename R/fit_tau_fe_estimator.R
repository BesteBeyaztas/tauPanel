fit_tau_fe_estimator <- function(
    y,
    x,
    id,
    N = 2000,
    kk = 2,
    tt = 20,
    rr = 2,
    approximate = 0,
    cluster_start = TRUE,
    m_cluster = NULL,
    seed = NULL
) {
  fit <- fit_fast_tau_scaled(
    x = x,
    y = y,
    N = N,
    kk = kk,
    tt = tt,
    rr = rr,
    approximate = approximate,
    seed = seed,
    id = id,
    cluster_start = cluster_start,
    m_cluster = m_cluster
  )
  
  list(beta = as.numeric(fit$beta), fit = fit)
}
