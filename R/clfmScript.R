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
#clfm(Triangle, selected)
#stop()
alpha <- tryCatch(clfmAlpha(Triangle, selected), # user can supply tolerance
                  error = function(e) e)
if (inherits(alpha, "error")) stop("Error solving for alpha")

CL <- chainladder(Triangle, delta = alpha)
clsoln <- CL

parms <- clParmsLogLinear(CL)
clfmsoln <- futureTriangles(Triangle, parms)

mcl <- MackChainLadder(Triangle)
print(clfmsoln)
print(mcl)

parms <- clParmsMack(CL)
clfmsoln <- futureTriangles(Triangle, parms)

mcl <- MackChainLadder(Triangle, est.sigma = "Mack")
print(clfmsoln)
print(mcl)

parms <- clParmsLinearCV(CL)
clfmsoln <- futureTriangles(Triangle, parms)

print(clfmsoln)
print(mcl)

parms <- clParms(CL)

parms <- clParms(CL, est.sigma = "Mack")
clfmsoln <- futureTriangles(Triangle, parms)

mcl <- MackChainLadder(Triangle, est.sigma = "Mack")
print(clfmsoln)
print(mcl)

parms <- clParms(CL, est.sigma = "log-linear")
clfmsoln <- futureTriangles(Triangle, parms)

mcl <- MackChainLadder(Triangle, est.sigma = "log-linear")
print(clfmsoln)
print(mcl)

#row.names(parms) <- c(colnames(Triangle)[-1L], "Inf")
#row.names(parms) <- colnames(Triangle)[-1L]
#parms[9, 2] <- .01
#parms[9, 3] <- 40
clparms <- parms
tail <- 1.0
tail.se <- 0.0
tail.sigma <- 0.0
tail.alpha = 1
clfmsoln <- futureTriangles(Triangle, parms)
summary(clfmsoln)
#summaryMack(clfmsoln)
print(clfmsoln)
mcl <- MackChainLadder(Triangle)
print(mcl)
sqrt(clfmsoln$FT.d^2 + clfmsoln$FT.g^2)
mcl$Mack.S.E
mcl$Total.Mack.S.E
sqrt(mcl$Total.ParameterRisk^2 + mcl$Total.ProcessRisk^2)
names(mcl)
round(clfmsoln$ft.d[,-11],0)
round(mcl$Mack.ParameterRisk, 0)
round(clfmsoln$ft.g[,-11],0)
round(mcl$Mack.ProcessRisk, 0)
round(clfmsoln$FT.d[-11], 0)
round(mcl$Total.ParameterRisk, 0)
round(clfmsoln$FT.g[-11], 0)
round(mcl$Total.ProcessRisk, 0)
mcl$Total.ProcessRisk
