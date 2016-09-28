#   futureTriangle.expd
mostCurrentDiagonal <- ChainLadder::getLatestCumulative
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
  
  ft <- ft.d <- ft.g <- matrix(NA, nrow = (m <- nrow(Triangle)), ncol = n1,
               dimnames = list(names(mcd), c( colnames(Triangle), "Inf")))

  
  # Project each observation/row
  for (i in 1:length(mcd)) {
    target <- targetnms[J0[i]]
    ft[i, target] <- mcd[i] * f[target]
    ft.d[i, target] <- f.se2[target] * mcd[i]^2
    ft.g[i, target] <- sigma2[target] * mcd[i]^alpha[target]
    if (J0[i] < n) {
      for (j in (J0[i]+1):n) {
        target <- targetnms[j]
        srce <- sourcenms[j]
        ft[i, target] <- ft[i, srce] * f[target]
        ft.d[i, target] <- ft[i, srce]^2 * f.se2[target] + f2[target] * ft.d[i, srce] +
          f.se2[target] * ft.d[i, srce]
        ft.g[i, target] <- #expectedValueXtoTheAlpha(ft[i, srce], alpha[target]) * sigma
          (ft[i, srce]^alpha[target] + .5*alpha[target]*(alpha[target]-1)*ft.g[i, srce]) *
          sigma2[target] + f2[target] * ft.g[i, srce]
      }
    }
  }
  # Project the total row
  FT <- FT.d <- FT.g <- structure(numeric(n1), names = colnames(ft))
  
  list(ft = ft, ft.d = sqrt(ft.d), ft.g = sqrt(ft.g), 
       FT = FT, FT.d = sqrt(FT.d), FT.g = sqrt(FT.g),
       clfmparms = data.frame(
    f = f, f.se = f.se, sigma = sigma, df = df, alpha = alpha, 
    source = sourcenms, target = targetnms))
}