#   futureTriangle.expd
mostCurrentDiagonal <- ChainLadder::getLatestCumulative
futureTriangles <- function(Triangle, clparms, 
                                 tail = 1.0, tail.se = 0, tail.sigma = 0, tail.alpha = 1) {
  f <- clparms$f
  f.se <- clparms$f.se
  sigma <- clparms$sigma
  alpha <- clparms$alpha
  df 	<- clparms$df
  n <- ncol(Triangle)
  if (length(f) < n) {
    f     <- c(f,     tail)
    f.se  <- c(f.se,  tail.se)
    sigma <- c(sigma, tail.sigma)
    alpha <- c(alpha, tail.alpha)
    df 	  <- c(df,    NA)
  }
  sigma2 <- sigma^2
  f2 <- f^2
  f.se2 <- f.se^2
	# All are named by the target node. "Inf" means last target is ultimate
  targetnms <- names(f)<-names(f.se)<-names(sigma)<-names(sigma2)<-
    names(f2)<-names(f.se2)<-names(alpha)<-c(row.names(clparms), "Inf")
	# If I'm assuming sourcenms = targetnms but for last, is that a problem?
  sourcenms <- colnames(Triangle)
  
  mcd <- getLatestCumulative(Triangle)
  J0 <- attr(mcd, "latestcol") # column in Triangle where each mcd resides

  n1 <- n+1
  
	# NA's will be in non-future time periods
  ft <- ft.d <- ft.g <- matrix(as.numeric(NA), nrow = (m <- nrow(Triangle)), ncol = n1,
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
          (ft[i, srce]^alpha[target] + 
             .5 * alpha[target] * (alpha[target]-1) * ft[i, srce]^(alpha[target] - 2) * ft.g[i, srce]) *
          sigma2[target] + f2[target] * ft.g[i, srce]
      }
    }
  }

  # Project stats for the total row
  FT.d <- FT.g <- structure(numeric(n1), names = colnames(ft))
  FT <- colSums(ft, na.rm = TRUE)
  FT[apply(ft, 2, function(x) all(is.na(x)))] <- NA
  # The total of the projected amounts in each column can be 
  # equivalently reproduced by recursively developing the sum of
  #   A. the previous total and B. the sum of the mcd's at that age.
  # The sum A+B is named TotalAmtDeveloped
	# First, calculate B; i.e., total up most current diagonal 
	#		at source of each development period
  mcdsourceage <- attr(mcd, "colnames")
  B <- sapply(sourcenms, function(s) sum(mcd[mcdsourceage == s]))
	# Then ..
  TotalAmtDeveloped <- cumsum(c(B, 0) + colSums(ft, na.rm = TRUE))
  
	# Using TotalAmtDeveloped, run through same formulas as above.
  for (i in 1:n) {
    srce <- sourcenms[i]
    target <- targetnms[i]
    FT.d[target] <- (TotalAmtDeveloped[i])^2 * f.se2[target] + f2[target] * FT.d[srce] +
      f.se2[target] * FT.d[srce]
    FT.g[target]  <- #expectedValueXtoTheAlpha(ft[i, srce], alpha[target]) * sigma
      (TotalAmtDeveloped[i]^alpha[target] + 
         .5 * alpha[target] * (alpha[target]-1) * 
          TotalAmtDeveloped[i]^(alpha[target] - 2) * FT.g[srce]) * sigma2[target] + 
          f2[target] * FT.g[srce]
  }
  
  list(ft = ft, ft.d = sqrt(ft.d), ft.g = sqrt(ft.g), 
       FT = FT, FT.d = sqrt(FT.d), FT.g = sqrt(FT.g),
       clfmparms = data.frame(f = f, f.se = f.se, sigma = sigma, 
														  df = df, alpha = alpha), 
			 source = sourcenms, target = targetnms)
}