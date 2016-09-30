#   futureTriangle.expd
mostCurrentDiagonal <- ChainLadder::getLatestCumulative
futureTriangles <- function(Triangle, clparms, 
                                 tail = 1.0, tail.se = 0, tail.sigma = 0, tail.alpha) {
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
    if (missing(tail.alpha)) tail.alpha <- tail(alpha, 1)
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
  TotalAmtDeveloped <- c(B, 0) + colSums(ft, na.rm = TRUE)
  
	# Using TotalAmtDeveloped, run through same formulas as above.
  for (i in 1:n) {
    srce <- sourcenms[i]
    target <- targetnms[i]
    FT.d[target] <- (TotalAmtDeveloped[srce])^2 * f.se2[target] + f2[target] * FT.d[srce] +
      f.se2[target] * FT.d[srce]
    FT.g[target]  <- #expectedValueXtoTheAlpha(ft[i, srce], alpha[target]) * sigma
      (TotalAmtDeveloped[srce]^alpha[target] + 
         .5 * alpha[target] * (alpha[target]-1) * 
          TotalAmtDeveloped[srce]^(alpha[target] - 2) * FT.g[srce]) * sigma2[target] + 
          f2[target] * FT.g[srce]
  }
  
  structure(
    list(ft = ft, ft.d = sqrt(ft.d), ft.g = sqrt(ft.g), 
         FT = FT, FT.d = sqrt(FT.d), FT.g = sqrt(FT.g),
         mcd = mcd,
         clfmparms = data.frame(f = f, f.se = f.se, sigma = sigma, 
				  										  df = df, alpha = alpha), 
         tail = tail, tail.se = tail.se, tail.sigma = tail.sigma, tail.alpha = tail.alpha,
			   source = sourcenms, target = targetnms),
    class = "clfm")
}
summary.clfm <- function(clfmsoln) {
  exh <- data.frame(
    origin = names(clfmsoln$mcd),
    latest = clfmsoln$mcd,
    ldf = cumprod(rev(clfmsoln$clfmparms$f)),
    ultimate = clfmsoln$ft[, "Inf"],
    se = sqrt(clfmsoln$ft.d[, "Inf"]^2 + clfmsoln$ft.g[, "Inf"]^2)
  )
  exh <- rbind(exh, data.frame(
    origin = "Sum",
    latest = sum(clfmsoln$mcd),
    ldf = as.numeric(NA),
    ultimate = sum(clfmsoln$ft[, "Inf"]),
    se = sqrt(clfmsoln$FT.d["Inf"]^2 + clfmsoln$FT.g["Inf"]^2)
    ))
  names(exh)[1L] <- attr(clfmsoln$mcd, "rowsname")
  exh
}
summaryMack <- function(x, ...) {UseMethod("summaryMack")}
summaryMack.MackChainLadder <- ChainLadder:::summary.MackChainLadder
summaryMack.clfm <- function(x, ...){
  
  ## Summarise my results
  Latest <- x$mcd
  
  ex.origin.period <- Latest!=0
  
  Ultimate <- x$ft[, "Inf"]
  Dev.To.Date <- Latest/Ultimate
  IBNR <- Ultimate-Latest
  Mack.S.E <- sqrt(x$ft.d[, "Inf"]^2 + x$ft.g[, "Inf"]^2)
  CV <- Mack.S.E/(Ultimate-Latest)
  
  ByOrigin <- data.frame(Latest, Dev.To.Date, Ultimate, IBNR, Mack.S.E, CV)
  names(ByOrigin)[6]="CV(IBNR)"
  ByOrigin <- ByOrigin[ex.origin.period,]
  
  Totals <-  c(sum(Latest,na.rm=TRUE),
               sum(Latest,na.rm=TRUE)/sum(Ultimate,na.rm=TRUE),
               sum(Ultimate,na.rm=TRUE),
               sum(IBNR,na.rm=TRUE), sqrt(x$FT.d["Inf"]^2 + x$FT.g["Inf"]^2),
               sqrt(x$FT.d["Inf"]^2 + x$FT.g["Inf"]^2)/sum(IBNR,na.rm=TRUE)
  )
  # Totals <- c(Totals, round(x[["Total.Mack.S.E"]]/sum(res$IBNR,na.rm=TRUE),2))
  Totals <- as.data.frame(Totals)
  
  colnames(Totals)=c("Totals")
  rownames(Totals) <- c("Latest:","Dev:","Ultimate:",
                        "IBNR:","Mack S.E.:",
                        "CV(IBNR):")
  
  output <- list(ByOrigin=ByOrigin, Totals=Totals)
  return(output)
}
print.clfm <- function(x, ...){
  
  ## Summarise my results
  Latest <- x$mcd
  
  ex.origin.period <- Latest!=0
  
  Ultimate <- x$ft[, "Inf"]
  Dev.To.Date <- Latest/Ultimate
  IBNR <- Ultimate-Latest
  Mack.S.E <- sqrt(x$ft.d[, "Inf"]^2 + x$ft.g[, "Inf"]^2)
  CV <- Mack.S.E/(Ultimate-Latest)
  
  ByOrigin <- data.frame(Latest, Dev.To.Date, Ultimate, IBNR, Mack.S.E, CV)
  names(ByOrigin)[6]="CV(IBNR)"
  ByOrigin <- ByOrigin[ex.origin.period,]
  print(format(ByOrigin, big.mark = ",", digits = 3),...)
  
  Totals <-  c(sum(Latest,na.rm=TRUE),
               sum(Latest,na.rm=TRUE)/sum(Ultimate,na.rm=TRUE),
               sum(Ultimate,na.rm=TRUE),
               sum(IBNR,na.rm=TRUE), sqrt(x$FT.d["Inf"]^2 + x$FT.g["Inf"]^2),
               sqrt(x$FT.d["Inf"]^2 + x$FT.g["Inf"]^2)/sum(IBNR,na.rm=TRUE)
  )
  # Totals <- c(Totals, round(x[["Total.Mack.S.E"]]/sum(res$IBNR,na.rm=TRUE),2))
  Totals <- as.data.frame(Totals)
  
  colnames(Totals)=c("Totals")
  rownames(Totals) <- c("Latest:","Dev:","Ultimate:",
                        "IBNR:","Mack S.E.:",
                        "CV(IBNR):")
  Totals[1:6,] <- formatC(Totals[1:6,], big.mark=",",digits=2,format="f")
  cat("\n")
  print(Totals, quote=FALSE)
#  
#  output <- list(ByOrigin=ByOrigin, Totals=Totals)
#  return(output)
}
