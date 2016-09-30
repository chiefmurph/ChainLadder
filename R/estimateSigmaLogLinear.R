# LogLinear estimate of ...
estimateSigmaLogLinear <- function(sigma) {
  isna <- is.na(sigma)
  if (sum(!isna) == 1) {
    stop(paste("Too few (1) link ratios for fitting 'loglinear' model to estimate sigma_n.\n",
               "Re-run with 'est.sigma = \"Mack\" or \"cvmethod\"'"))
  }
  ## estimate sigma[n-1] via log-linear regression
  sig.model <- suppressWarnings(ChainLadder:::estimate.sigma(sigma))
  p.value.of.model <- tryCatch(summary(sig.model$model)$coefficient[2,4],
                               error = function(e) e)
  if (inherits(p.value.of.model, "error"))
    stop("p-value of 'loglinear' sigma model indeterminable.")
  if(p.value.of.model > 0.05)
        warning(paste0("'loglinear' model to estimate sigma may not be appropriate. p-value = ",
                       p.value.of.model))
  sig.model$sigma
}
estimateSigmaMack <- function(sigma) {
  for(i in which(is.na(sigma))) {   # usually i = n - 1
    sigma[i] <- sqrt(abs(min((sigma[i - 1]^4/sigma[i - 2]^2),
                             min(sigma[i - 2]^2, sigma[i - 1]^2))))
  }
  sigma
}
estimateSigmaLinearCV <- function(sigma, f) {
  w <- which(is.na(sigma))
  x <- f - 1
  m <- lm(sigma ~ x)
  sigma[w] <- predict(m, newdata = data.frame(x = x[w]))
  sigma
}
