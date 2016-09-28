clParms <- function(clSoln) {
  smmry <- suppressWarnings(lapply(clSoln$Models, summary))
  f <- sapply(smmry, function(x) x$coef["x","Estimate"])
  f.se <- sapply(smmry, function(x) x$coef["x","Std. Error"])
  sigma <- sapply(smmry, function(x) x$sigma)
  df <- sapply(smmry, function(x) x$df[2L])
  alpha <- clSoln$delta
  # Check for NA/NaN variance parameters
  w <- which(is.na(sigma))
  if (w) {
    v <- mean(sigma[-w] / (f[-w] - 1))
    sigma[w] <- v * (f[w] - 1)
    f.se[w] <- sapply(w, function(s){
      dat <- clSoln$Models[[s]]$model
      sigma[s] / sqrt(sum(dat$x^(2-alpha[s])))
    })
  }
  data.frame(f, f.se, sigma, df, alpha)
}
fminus1.cv <- function(macksoln) {
  x <- mackParms(macksoln)
  x$f.se / (x$f - 1)
}
