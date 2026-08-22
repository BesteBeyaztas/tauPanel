FastTau <- function(x, y, N=500, kk=2, tt=5, rr=2, approximate=0, seed=NULL,
                    id = NULL, cluster_start = FALSE, m_cluster = NULL)
{
  if (tt<1) stop("parameter tt should be at least 1")
  
  had_random_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  
  if (had_random_seed) {
    old_random_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }
  
  on.exit(
    {
      if (had_random_seed) {
        
        assign(
          ".Random.seed",
          old_random_seed,
          envir = .GlobalEnv
        )
        
      } else if (
        exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      ) {
        
        rm(
          ".Random.seed",
          envir = .GlobalEnv
        )
      }
    },
    add = TRUE
  )
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  x <- as.matrix(x)
  y <- as.vector(y)
  
  n <- nrow(x)
  p <- ncol(x)
  
  if (!is.null(id)) {
    id <- as.factor(id)
  }
  
  if (cluster_start) {
    if (is.null(id)) {
      stop("If cluster_start=TRUE, an id vector must be supplied.")
    }
    
    id_levels <- levels(id)
    G <- length(id_levels)
    
    if (G < p) {
      stop(
        "The number of individuals must be at least the number ",
        "of regression coefficients for cluster-based starts."
      )
    }
    
    if (is.null(m_cluster)) {
      m_cluster <- p
    }
    
    if (length(m_cluster) != 1 ||
        !is.finite(m_cluster) ||
        abs(m_cluster - round(m_cluster)) > sqrt(.Machine$double.eps) ||
        m_cluster < p ||
        m_cluster > G) {
      stop(
        "m_cluster must be an integer between the number of ",
        "regression coefficients and the number of individuals."
      )
    }
    
    m_cluster <- as.integer(m_cluster)
  }
  
  draw_start_rows <- function() {
    if (!cluster_start) {
      return(randomset(n, p))
    }
    
    chosen_ids <- sample(
      id_levels,
      size = m_cluster,
      replace = FALSE
    )
    
    selected_rows <- vapply(
      chosen_ids,
      function(current_id) {
        available_rows <- which(
          as.character(id) == current_id
        )
        sample(available_rows, size = 1)
      },
      integer(1)
    )
    
    unname(selected_rows)
  }
  
  c1 <- .4046
  b1 <- .5
  c2 <- 1.09
  b2 <- .1278
  
  RWLStol <- 1e-11
  
  bestbetas <- matrix(0, p, tt)
  bestscales <- 1e20 * rep(1, tt)
  besttauscales <- 1e20 * rep(1, tt)
  worsti <- 1
  rworst <- y
  
  for (i in 1:N) {
    singular <- 1; itertest <- 1
    while (singular==1 && itertest<100) {
      ranset <- draw_start_rows()
      xj <- x[ranset, , drop = FALSE]
      yj <- y[ranset]
      
      qx <- qr(xj)
      
      if (qx$rank < p) {
        singular <- 1
      } else {
        bj <- as.matrix(qr.coef(qx, yj))
        singular <- any(!is.finite(bj))
      }
      itertest <- itertest + 1
    }
    if (singular == 1) {
      stop("Too many degenerate subsamples")
    }
    
    if (kk > 0) {
      tmp <- IWLSiteration(x, y, bj, 0, kk, RWLStol, b1, c1, c2)
      betarw <- tmp$betarw
      resrw <- y - x %*% betarw
      scalerw <- tmp$scalerw
    }
    else {
      betarw <- bj
      resrw <- y - x %*% betarw
      scalerw <- median(abs(resrw))/.6745
    }
    
    if (i > 1) LTMvec = LTMvec + abs(resrw)
    else LTMvec = abs(resrw)
    
    
    if (!approximate) {
      scaletest1 <- mean(rhoOpt(resrw / bestscales[worsti],c1)) < b1
      scaletest2 <- sum(rhoOpt(resrw / bestscales[worsti],c2)) < sum(rhoOpt(rworst/bestscales[worsti],c2))
      if (scaletest1 || scaletest2) {
        snew <- Mscale(resrw, b1, c1, scalerw)
        taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
        if (taunew < besttauscales[worsti]) {
          besttauscales[worsti] <- taunew
          bestscales[worsti] <- snew
          bestbetas[,worsti] <- betarw
          worsti <- which.max(besttauscales) 
          rworst <- y - x %*% bestbetas[,worsti]
        }
      }
    }
    else {
      snew = scalerw;
      if (rr>0) {
        for (kstep in 1:rr) { 
          snew <- sqrt( snew^2 * mean( rhoOpt(resrw/snew,c1) ) / b1 )
        }
      }
      taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
      if (taunew < besttauscales[worsti]) {
        besttauscales[worsti] <- taunew
        bestscales[worsti] <- snew
        bestbetas[,worsti] <- betarw
        worsti <- which.max(besttauscales) 
        rworst <- y - x %*% bestbetas[,worsti]
      }
    }
  }
  
  if (!cluster_start) {
    
    IXLTM <- order(LTMvec, decreasing = TRUE)
    singular <- 1 
    extrasize <- p
    
    while (singular == 1) {
      xs <- x[IXLTM[1:extrasize], , drop = FALSE]
      ys <- y[IXLTM[1:extrasize]]
      qx <- qr(xs)
      
      if (qx$rank < p) {
        singular <- 1
      } else {
        bbeta <- as.matrix(qr.coef(qx, ys))
        singular <- any(!is.finite(bbeta))
      }
      
      extrasize <- extrasize + 1
    }
    
  } else {
    
    LTM_id <- tapply(as.vector(LTMvec), id, mean, na.rm = TRUE)
    bad_ids <- names(sort(LTM_id, decreasing = TRUE))
    
    singular <- 1
    extrasize <- m_cluster
    
    while (singular == 1 && extrasize <= length(bad_ids)) {
      chosen_ids <- bad_ids[seq_len(extrasize)]
      
      rows_extra <- vapply(
        chosen_ids,
        function(current_id) {
          current_rows <- which(as.character(id) == current_id)
          current_rows[which.max(LTMvec[current_rows])]
        },
        integer(1)
      )
      
      xs <- x[rows_extra, , drop = FALSE]
      ys <- y[rows_extra]
      qx <- qr(xs)
      
      if (qx$rank < p) {
        singular <- 1
      } else {
        bbeta <- as.matrix(qr.coef(qx, ys))
        singular <- any(!is.finite(bbeta))
      }
      
      extrasize <- extrasize + 1
    }
    
    if (singular == 1) {
      bbeta <- as.matrix(qr.coef(qr(x), y))
    }
  }
  
  if (kk > 0) {
    tmp <- IWLSiteration(x, y, bbeta, 0, kk, RWLStol, b1, c1, c2)
    betarw <- tmp$betarw
    resrw <- y - x %*% betarw
    scalerw <- tmp$scalerw
  }
  else {
    betarw <- bbeta
    resrw <- y - x %*% betarw
    scalerw <- median(abs(resrw))/.6745
  }
  
  if (!approximate) {
    scaletest1 <- mean(rhoOpt(resrw / bestscales[worsti],c1)) < b1
    scaletest2 <- sum(rhoOpt(resrw / bestscales[worsti],c2)) < sum(rhoOpt(rworst/bestscales[worsti],c2))
    if (scaletest1 || scaletest2) {
      snew <- Mscale(resrw, b1, c1, scalerw)
      taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
      if (taunew < besttauscales[worsti]) {
        besttauscales[worsti] <- taunew
        bestscales[worsti] <- snew
        bestbetas[,worsti] <- betarw
        worsti <- which.max(besttauscales) 
        rworst <- y - x %*% bestbetas[,worsti]
      }
    }
  }
  else {
    snew = scalerw;
    if (rr>0) {
      for (kstep in 1:rr) { 
        snew <- sqrt( snew^2 * mean( rhoOpt(resrw/snew,c1) ) / b1 )
      }
    }  
    taunew <- snew * sqrt(mean(rhoOpt(resrw/snew,c2)))
    if (taunew < besttauscales[worsti]) {
      besttauscales[worsti] <- taunew
      bestscales[worsti] <- snew
      bestbetas[,worsti] <- betarw
      worsti <- which.max(besttauscales) 
      rworst <- y - x %*% bestbetas[,worsti]
    }
  }
  
  superbesttauscale <- 1e20
  
  for (i in 1:tt) {
    tmp <- IWLSiteration(x, y, bestbetas[,i], bestscales[i], 500, RWLStol, b1, c1, c2)
    resrw <- y - x %*% tmp$betarw
    tauscalerw <- tmp$scalerw * sqrt(mean(rhoOpt(resrw/tmp$scalerw,c2)))
    if (tauscalerw < superbesttauscale) {
      superbesttauscale <- tauscalerw
      superbestbeta <- tmp$betarw
      superbestscale <- tmp$scalerw
    }
  }
  
  superbestscale <- Mscale(y - x%*%superbestbeta, b1, c1, superbestscale)
  superbesttauscale <- superbestscale * sqrt(mean(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2)))
  
  betaLS <- as.matrix(qr.coef(qr(x),y))
  resLS <- y - x %*% betaLS
  scaleLS <- median(abs(resLS))/.6745  
  scaletest1 <- mean(rhoOpt(resLS / superbestscale,c1)) < b1
  scaletest2 <- sum(rhoOpt(resLS / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
  if (scaletest1 || scaletest2) {
    snew <- Mscale(resLS, b1, c1, scaleLS)
    taunew <- snew * sqrt(mean(rhoOpt(resLS/snew,c2)))
    if (taunew < superbesttauscale) {
      superbestscale <- snew
      superbestbeta <- betaLS
      superbesttauscale <- taunew
    }   
  }
  
  IXmed <- order(abs(y - median(y)))
  xhalf <- x[IXmed[1:floor(n/2)],]
  yhalf <- y[IXmed[1:floor(n/2)]]
  bbeta <- as.matrix(qr.coef(qr(xhalf),yhalf))
  tmp <- IWLSiteration(x, y, bbeta, 0, 10, RWLStol, b1, c1, c2)
  betaHB <- tmp$betarw
  resHB <- y - x %*% betaHB
  scaleHB <- tmp$scalerw
  scaletest1 <- mean(rhoOpt(resHB / superbestscale,c1)) < b1
  scaletest2 <- sum(rhoOpt(resHB / superbestscale,c2)) < sum(rhoOpt((y - x%*%superbestbeta)/superbestscale,c2))
  if (scaletest1 || scaletest2) {
    snew <- Mscale(resHB, b1, c1, scaleHB)
    taunew <- snew * sqrt(mean(rhoOpt(resHB/snew,c2)))
    if (taunew < superbesttauscale) {
      superbestbeta <- betaHB
      superbesttauscale <- taunew
    }   
  }
  
  
  return(list( beta = superbestbeta, scale = superbesttauscale/sqrt(b2) ))
  
}
