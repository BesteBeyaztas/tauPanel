randomset <- function(tot,nel) {
  ranset <- rep(0,nel)
  for (j in 1:nel) {
    num <- ceiling(runif(1)*tot)
    if (j > 1) {
      while (any(ranset==num)) 
        num <- ceiling(runif(1)*tot)
    }
    ranset[j] <- num
  }
  return(ranset)
}
