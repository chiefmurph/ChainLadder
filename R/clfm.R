# clfm.r
# Calculate the mse per the CLFM paper, Variance vol 6 issue 2
#   with the modification of using Delta method r.t. psi function
#   when calculating process risk recursively.
source(file.path("R", "clParms.R"))
source(file.path("R", "futureTriangles.R"))
clfmAlpha <- function(Triangle, ...) {
  x <- CLFMdelta(Triangle, ...)
  names(x) <- colnames(Triangle)[-1L]
  x
}

clfm <- function(Triangle, selected, ...) {
  
  Triangle <- ChainLadder:::checkTriangle(Triangle)
  
  alpha <- tryCatch(clfmAlpha(Triangle, selected, ...), # user can supply tolerance
                    error = function(e) e)
  if (inherits(alpha, "error")) stop("Error solving for alpha")
  
  ##################
  # Create a list of all the models
  CL <- chainladder(Triangle, delta = alpha)
  parms <- clParms(CL)
   
  futureTriangles(Triangle, parms)
}
otherstuff <- function() {
  FullTriangle <- predict.ChainLadder(list(Models=CL[["Models"]], Triangle = Triangle))
#  f <- rep(1, n - 1)
#  f.se <- rep(0, n - 1)
#  sigma <- rep(0, n - 1)
  
  ##################
  smmry <- lapply(CL[["Models"]], summary)
  f <- sapply(smmry, function(x) x$coef["x", "Estimate"])
  f.se <- sapply(smmry, function(x) x$coef["x", "Std. Error"])
  sigma <- sapply(smmry, function(x) x$sigma)
  df <- sapply(smmry, function(x) x$df[2L])
  
  ##################
  tolerance <- .Machine$double.eps
  perfect.fit <- (df > 1) & (f.se < tolerance)
  w <- which(perfect.fit)
  if (length(w)) {
    warn <- "Information: essentially no variation in development data for period(s):\n"
    nms <- colnames(FullTriangle)
    periods <- paste0("'", paste(nms[w], nms[w+1], sep = "-"), "'")
    warn <- c(warn, paste(periods, collapse = ", "))
    cat(warn, "\n")
  }
  
  ##################
  for(i in which(isna)){   # usually i = n - 1
    sigma[i] <- sqrt(abs(min((sigma[i - 1]^4/sigma[i - 2]^2),
                             min(sigma[i - 2]^2, sigma[i - 1]^2))))
    f.se[i] <- sigma[i] / sqrt(weights[1,i] * FullTriangle[1,i]^alpha[i])
  }

  ##################  
  FullTriangle.paramrisk <- matrix(0, m, n)
  for (k in 1:(n-1)) {
    for (i in (m-k+1):m) {
      FullTriangle.procrisk[i, k+1] <- 0
    }
  }
  
  ##################
  FullTriangle.procrisk <- matrix(0, m, n)
  FullTriangle.procrisk <- fullTriangleProcRisk(macksoln)
}

LM.MackHeuristic.ATA <- function(LM,k){
  if (k<3) 
    stop("Not enough prior periods to implement Mack's heuristic at age",k)
  sigma <- sqrt(min(LM[[k-1]]$sigma^4/LM[[k-2]]$sigma^2,
                    min(LM[[k-2]]$sigma^2, LM[[k-1]]$sigma^2))
  )    
  # That implements Mack's idea for sigma. The formula for sef is:
  #     sef^2 = sigma^2 / denominator
  # where denominator = sum of beginning losses ^ (2-alpha)
  #       in the calculation of this particular ata.
  # Those beginning loss values are the values in LM[[.]]$model in the
  #   column whose name is stored in LM[[.]]$xname
  sef <- sigma/sqrt(sum(LM[[k]]$model[LM[[k]]$xname]^(2-LM[[k]]$alpha)))
  # For alpha and df -- not aspects of Mack's original method --
  #   we will use linear interpolation on the prior two periods, 
  #   where the proportionality factor is based on the sigma's.
  lambda <- (sigma-LM[[k-2]]$sigma)/(LM[[k-1]]$sigma-LM[[k-2]]$sigma)
  df <- LM[[k-2]]$df + lambda * (LM[[k-1]]$df-LM[[k-2]]$df)
  # changed alpha from interpolation to prior alpha 3/26/10
  alpha <- LM[[k-1]]$alpha
  #    alpha <- LM[[k-2]]$alpha + lambda * (LM[[k-1]]$alpha-LM[[k-2]]$alpha)
  #    ffmInfoMsg(paste("Implementing Mack ata error heuristics at age ",
  # changed from ffmInfoMsg to trowarning because the order of the message
  # showing up on the console was reverse of actual order
  trowarning(paste("Implementing Mack ata error heuristics at age ",
                   LM[[k]]$age, ", alpha = prior alpha = ", alpha, sep=""))
  list(sef=sef,
       sigma=sigma,
       df=df,
       alpha=alpha)
}

