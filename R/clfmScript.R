# clfmScript
library(ChainLadder)
source(file.path("R", "clfm.R"))
#source(file.path("R", "clParms.R"))
#source(file.path("R", "futureTriangles.R"))
#clfmAlpha <- function(Triangle, ...) {
#  x <- CLFMdelta(Triangle, ...)
#  names(x) <- colnames(Triangle)[-1L]
#  x
#}
Triangle <- GenIns
#Triangle <- ChainLadder:::checkTriangle(Triangle)
selected <- attr(ata(Triangle), "vwtd")
clfm(Triangle, selected)
alpha <- tryCatch(clfmAlpha(Triangle, selected), # user can supply tolerance
                  error = function(e) e)
if (inherits(alpha, "error")) stop("Error solving for alpha")

CL <- chainladder(Triangle, delta = alpha)
parms <- clParms(CL)
#row.names(parms) <- c(colnames(Triangle)[-1L], "Inf")
#row.names(parms) <- colnames(Triangle)[-1L]
#parms[9, 2] <- .01
#parms[9, 3] <- 40
#clparms <- parms
#tail <- 1.0
#tail.se <- 0.0
#tail.sigma <- 0.0
future <- futureTriangles(Triangle, parms)
