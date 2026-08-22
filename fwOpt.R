fwOpt <- function(x, cc)
{
  tmp <- (-1.944 / cc^2 + 1.728 * x^2 / cc^4 - 0.312 * x^4 / cc^6 + 0.016 * x^6 / cc^8) / 3.25
  tmp[abs(x) < 2*cc] <- 1 / (3.25*cc^2)
  tmp[abs(x) > 3*cc] <- 0
  tmp
  
}