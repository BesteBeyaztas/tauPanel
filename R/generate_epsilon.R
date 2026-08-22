generate_epsilon <- function(N_, T_, distribution, sigma2_e, chi_df = 1) {
  
  if (distribution == "normal") {
    
    epsilon <- matrix(
      rnorm(N_ * T_, sd = sqrt(sigma2_e)),
      nrow = N_
    )
    
  } else if (distribution == "t") {
    
    epsilon <- matrix(rt(N_ * T_, df = 5), nrow = N_)
    
  } else if (distribution == "chi") {
    
    epsilon <- matrix(rchisq(N_ * T_, df = chi_df), nrow = N_)
    
  } else {
    
    stop("Unknown error distribution.")
  }
  
  return(epsilon)
}
