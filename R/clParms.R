clParms <- function(clsoln, est.sigma = "linearcv") {
  est.sigma <- match.arg(est.sigma, c(
    "linearcv", "Mack", "log-linear"
  ))
  smmry <- suppressWarnings(lapply(clsoln$Models, summary))
  f <- sapply(smmry, function(x) x$coef["x","Estimate"])
  f.se <- sapply(smmry, function(x) x$coef["x","Std. Error"])
  sigma <- sapply(smmry, function(x) x$sigma)
  # Check for NA/NaN variance parameters
  w <- which(is.na(sigma))
  if (w) {
    sigma <- switch(est.sigma,
                   linearcv = estimateSigmaLinearCV(sigma, f),
                   Mack = estimateSigmaMack(sigma),
                   `log-linear` = estimateSigmaLogLinear(sigma))
    f.se[w] <- sapply(w, function(s){
      dat <- clsoln$Models[[s]]$model
      sigma[s] / sqrt(sum(dat$x^(2-alpha[s])))
    })
  }
  df <- sapply(smmry, function(x) x$df[2L])
  alpha <- clsoln$delta
  data.frame(f, f.se, sigma, df, alpha)
}
fminus1.cv <- function(macksoln) {
  x <- mackParms(macksoln)
  x$f.se / (x$f - 1)
}
clParmsLogLinear <- function(clsoln) {
  smmry <- suppressWarnings(lapply(clsoln$Models, summary))
  f <- sapply(smmry, function(x) x$coef["x","Estimate"])
  f.se <- sapply(smmry, function(x) x$coef["x","Std. Error"])
  sigma <- sapply(smmry, function(x) x$sigma)
  # Check for NA/NaN variance parameter
  w <- which(is.na(sigma))
  if (w) {
    sigma <- estimateSigmaLogLinear(sigma)
    f.se[w] <- sapply(w, function(s){
      dat <- clsoln$Models[[s]]$model
      sigma[s] / sqrt(sum(dat$x^(2-alpha[s])))
    })
  }
  df <- sapply(smmry, function(x) x$df[2L])
  alpha <- clsoln$delta
  data.frame(f, f.se, sigma, df, alpha)
}
clParmsMack <- function(clsoln) {
  smmry <- suppressWarnings(lapply(clsoln$Models, summary))
  f <- sapply(smmry, function(x) x$coef["x","Estimate"])
  f.se <- sapply(smmry, function(x) x$coef["x","Std. Error"])
  sigma <- sapply(smmry, function(x) x$sigma)
  # Check for NA/NaN variance parameter
  w <- which(is.na(sigma))
  if (w) {
    sigma <- estimateSigmaMack(sigma)
    f.se[w] <- sapply(w, function(s) {
      dat <- clsoln$Models[[s]]$model
      sigma[s] / sqrt(sum(dat$x^(2-alpha[s])))
    })
  }
  df <- sapply(smmry, function(x) x$df[2L])
  alpha <- clsoln$delta
  data.frame(f, f.se, sigma, df, alpha)
}
clParmsLinearCV <- function(clsoln) {
  smmry <- suppressWarnings(lapply(clsoln$Models, summary))
  f <- sapply(smmry, function(x) x$coef["x","Estimate"])
  f.se <- sapply(smmry, function(x) x$coef["x","Std. Error"])
  sigma <- sapply(smmry, function(x) x$sigma)
  # Check for NA/NaN variance parameter
  w <- which(is.na(sigma))
  if (w) {
    sigma <- estimateSigmaLinearCV(sigma, f)
    f.se[w] <- sapply(w, function(s) {
      dat <- clsoln$Models[[s]]$model
      sigma[s] / sqrt(sum(dat$x^(2-alpha[s])))
    })
  }
  df <- sapply(smmry, function(x) x$df[2L])
  alpha <- clsoln$delta
  data.frame(f, f.se, sigma, df, alpha)
}
