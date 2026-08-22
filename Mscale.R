Mscale <- function(u, b, c, initialsc) 
{
  if (initialsc == 0) {
    initialsc <- median(abs(u)) / .6745
  }
  
  initialsc <- max(initialsc, 1e-8)
  maxit <- 100
  sc <- initialsc
  i <- 0 
  eps <- 1e-10
  err <- 1
  while  (( i < maxit ) & (err > eps)) {
    sc2 <- sqrt( sc^2 * mean(rhoOpt(u/sc,c)) / b)
    err <- abs(sc2/sc - 1)
    sc <- sc2
    i <- i+1
  }
  
  return(sc)
  
}