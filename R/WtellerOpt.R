WtellerOpt <- function(x, cc)
{
  tmp <- (3.584 - 0.864 * x^4 / cc^4 + 0.208 * x^6 / cc^6 - 0.012 * x^8 / cc^8) / 3.25
  tmp[abs(x) < 2*cc] <- 0
  tmp[abs(x) > 3*cc] <- 2
  tmp
  
}
