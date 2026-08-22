IWLSiteration <- function(x, y, inib, iniscale, maxiter, tol, b, c1, c2)
{  
  n <- nrow(x)
  p <- ncol(x)
  
  res <- y - x %*% inib
  if (iniscale == 0) {
    scale <- median(abs(res)) / .6745
  } else {
    scale <- iniscale
  }
  
  scale <- max(scale, 1e-8)
  
  oldbeta <- inib
  newbeta <- inib
  
  betadiff <- 2*tol
  iter <- 0
  while ((betadiff > tol) && (iter < maxiter)) {
    scale <- sqrt(scale^2 * mean(rhoOpt(res / scale, c1)) / b)
    
    if (!is.finite(scale) || scale <= 0) {
      scale <- 1e-8
    }
    scaledres <- res/scale
    Wn.teller <- sum(WtellerOpt(scaledres,c2))
    Wn.noemer <- sum(psixOpt(scaledres,c1))
    Wn <- Wn.teller / Wn.noemer
    weights <- (Wn * fwOpt(scaledres,c1) + fwOpt(scaledres,c2))
    sqweights <- sqrt(weights)
    xw <- x * as.vector(sqweights)
    yw <- y * sqweights
    newbeta <- qr.coef(qr(xw),yw)
    if (any(!is.finite(newbeta))) {
      newbeta <- inib
      scale <- iniscale
      break
    }
    betadiff <- sqrt(sum((oldbeta - newbeta)^2))
    res <- y - x %*% newbeta
    oldbeta <- newbeta
    iter <- iter + 1
  }
  
  return( list( betarw = newbeta, scalerw = scale ) )
  
}
