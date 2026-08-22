fast_tau_panel <- function(Y, X, effect = c("fe", "re"), ...) {
  effect <- match.arg(effect)
  
  if (effect == "fe") {
    fast_tau_fe(Y = Y, X = X, ...)
  } else {
    fast_tau_re(Y = Y, X = X, ...)
  }
}
