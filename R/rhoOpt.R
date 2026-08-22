rhoOpt <- function(x, cc)
{
  tmp <- x^2 / 2 / (3.25*cc^2)
  tmp2 <- (1.792 - 0.972 * x^2 / cc^2 + 0.432 * x^4 / cc^4 - 0.052 * x^6 / cc^6 + 0.002 * x^8 / cc^8) / 3.25
  tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
  tmp[abs(x) > 3*cc] <- 1
  tmp
  
}
