#   futureTriangle.expd
mostCurrentDiagonal <- ChainLadder::getLatestCumulative
futureTriangle.expd <- function(Triangle, clparms, tail = 1.0) {
  f <- clparms$f
  if (length(f) < (n<-ncol(Triangle))) f <- c(f, tail)
  ndx <- 1:n
  targetnms <- c(colnames(Triangle)[-1], "Inf")
  names(f) <- nms
  names(ndx) <- nms
  
  mcd <- getLatestCumulative(Triangle)
  ages <- nms[attr(mcd, "latestcol")] # name of the target node

  L <- lapply(1:length(mcd), function(i){
    cumprod(f[ndx[ages[i]]:n]) * mcd[i]
  })
  FT <- matrix(NA, nrow = (m <- nrow(Triangle)), ncol = n+1,
               dimnames = list(names(mcd),
                               c(colnames(Triangle)[1], names(f))))
  for (i in 1:m) {
    FT[i, names(L[[i]])] <- L[[i]]
  }
  FT
}

futureTriangle.param <- function(Triangle, clparms, 
                                 tail = 1.0, tail.se = 0, tail.sigma = 0) {
  f <- clparms$f
  f.se <- clparms$f.se
  sigma <- clparms$sigma
  if (length(f) < (n<-ncol(Triangle))) {
    f <- c(f, tail)
    f.se <- c(f, tail.se)
    sigma <- c(f, tail.sigma)
  }
  names(f)<-names(f.se)<-names(sigma) <- c(row.names(clparms), "Inf")
  
  n1 <- n+1
  
  ndx <- 1:n
  nms1 <- c(colnames(Triangle), "Inf")
  nms <- nms1[-1]
  names(f) <- nms
  names(ndx) <- nms
  
  mcd <- getLatestCumulative(Triangle)
  J0 <- attr(mcd, "latestcol")
  ages <- nms[attr(mcd, "latestcol")] # name of the target node
  age <- nms[attr(mcd, "latestcol")] # name of the target node
  
  FT <- matrix(NA, nrow = (m <- nrow(Triangle)), ncol = n+1,
               dimnames = list(names(mcd),
                               nms1))
  i = 1
  f.se[age[i]]
  FT[i, age[i]] <- 
  nms[i]
  nms[age[i]]
    for (i in 1:m) {
    jnms <- 
      FT[i, names(L[[i]])] <- L[[i]]
  }
  L <- lapply(1:length(mcd), function(i) {
    j <- ndx[ages[i]]:n
    x <- numeric(length(j))
    x[1L] <- f.se[j[1L]]^2 * mcd[i]^2
    for (t in j[-1L]) x[t] <- f.se[j[1L]]^2 * mcd[i]^2
  })
  FT
}

futureTriangles <- function(Triangle, clparms, 
                                 tail = 1.0, tail.se = 0, tail.sigma = 0, tail.alpha = 1) {
  f <- clparms$f
  f.se <- clparms$f.se
  sigma <- clparms$sigma
  alpha <- clparms$alpha
  n <- ncol(Triangle)
  if (length(f) < n) {
    f <- c(f, tail)
    f.se <- c(f.se, tail.se)
    sigma <- c(sigma, tail.sigma)
    alpha <- c(alpha, tail.alpha)
    df <- c(clparms$df, NA)
  }
  sigma2 <- sigma^2
  f2 <- f^2
  f.se2 <- f.se^2
  targetnms <- names(f)<-names(f.se)<-names(sigma)<-names(sigma2)<-
    names(f2)<-names(f.se2)<-names(alpha)<-c(row.names(clparms), "Inf")
  sourcenms <- colnames(Triangle)
  
  mcd <- getLatestCumulative(Triangle)
  J0 <- attr(mcd, "latestcol")

  n1 <- n+1
  
  FT <- FT.d <- FT.g <- matrix(NA, nrow = (m <- nrow(Triangle)), ncol = n1,
               dimnames = list(names(mcd), c( colnames(Triangle), "Inf")))

  
  # Project each observation/row
  for (i in 1:length(mcd)) {
    target <- targetnms[J0[i]]
    FT[i, target] <- mcd[i] * f[target]
    FT.d[i, target] <- f.se2[target] * mcd[i]^2
    FT.g[i, target] <- sigma2[target] * mcd[i]^alpha[target]
    if (J0[i] < n) {
      for (j in (J0[i]+1):n) {
        target <- targetnms[j]
        srce <- sourcenms[j]
        FT[i, target] <- FT[i, srce] * f[target]
        FT.d[i, target] <- FT[i, srce]^2 * f.se2[target] + f2[target] * FT.d[i, srce] +
          f.se2[target] * FT.d[i, srce]
        FT.g[i, target] <- #expectedValueXtoTheAlpha(FT[i, srce], alpha[target]) * sigma
          (FT[i, srce]^alpha[target] + .5*alpha[target]*(alpha[target]-1)*FT.g[i, srce]) *
          sigma2[target] + f2[target] * FT.g[i, srce]
      }
    }
  }
  list(FT = FT, FT.d = sqrt(FT.d), FT.g = sqrt(FT.g), clfmparms = data.frame(
    f = f, f.se = f.se, sigma = sigma, df = df, alpha = alpha, 
    source = sourcenms, target = targetnms))
}