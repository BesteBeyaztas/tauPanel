fit_tau_re_estimator <- function(
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
  components <- estimate_tau_re_components(
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
  
  tau_re_data <- panel_tau_re_quasi_demean(
    Y = Y,
    X = X,
    theta = components$theta,
    residual_matrix = components$residual_matrix
  )
  
  final_fit <- fit_fast_tau_scaled(
    x = tau_re_data$x,
    y = tau_re_data$y,
    N = N,
    kk = kk,
    tt = tt,
    rr = rr,
    approximate = approximate,
    seed = seed,
    id = tau_re_data$id,
    cluster_start = cluster_start,
    m_cluster = m_cluster
  )
  
  list(
    beta = as.numeric(final_fit$beta),
    sigma2_a = components$sigma2_a,
    sigma2_e = components$sigma2_e,
    theta = components$theta,
    fit = list(
      final = final_fit,
      preliminary = components$preliminary_fit,
      beta_initial = components$beta_initial,
      transformed_data = tau_re_data
    )
  )
}