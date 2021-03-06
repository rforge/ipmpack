\name{getIPMoutput}
\alias{getIPMoutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gets IPM output from a list of T matrices (only survival and size information) or full IPMs (T matrices + F matrices; the latter include sexual reproduction information).
}
\description{
Gets synthetic values including life expectancy, passage time, and, if a
fecundity matrix is available, population growth rate (lambda), stable stage distribution, reproductive output, etc. 
}
\usage{
getIPMoutput(TmatrixList, targetSize = c(), FmatrixList = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TmatrixList}{
    List of survival-size (T) matrices for which summary statistics desired.
}
  \item{targetSize}{
    Target size for passage time. If none is provided defaults to the median of the IPM meshpoints.
}
  \item{FmatrixList}{
    List of fecundity (F) matrices corresponding to the list of T matrices, facultative. When this information is included, population growth rate (lambda), and stable stage distribution will also be provided.
}
}
\value{
  \item{LE}{matrix of life expectancies, columns correspond to
  meshpoints, rows corresponding to each element of the list of T matrices}
  \item{pTime}{matrix of passage times to the targetSize from each of
    the meshpoints (columns) and for each element in the T matrix list
    (columns).}
  \item{lambda}{vector of population growth rates corresponding to value obtained
    combining each element of the list of T matrices with the
    corresponding element in the list of F matrices; if no F matrix list
    is provided, it returns a vector of NAs.}
  \item{stableStage}{matrix of stable stage distributions rows corresponding to values obtained.
    combining each element of the list of Tmatrices with the
    corresponding element in the list of Fmatrices; if no Fmatrix list
    is provided, this is a matrix of NAs.}
}
\references{
Zuidema, Jongejans, Chien, During & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{getIPMoutputDirect}} 
}
\examples{
dff <- generateData()

# Use package MCMCglmm to get Bayesian posteriors 
# NOT RUN ##
#svlist <- makePostSurvivalObjs(dataf = dff, explanatoryVariables = "size",
#meanB = rep(0, 3), varB = rep(1e+10), nitt = 5000)
#grlist <- makePostGrowthObjs(dataf = dff,  nitt = 5000)
#fvlist <- makePostFecObjs(dataf = dff, nitt = 5000)

# Use only first 10 of each list, for speed:
#Tlist <- makeListTmatrix(growObjList = grlist[1:10], survObjList = svlist[1:10], 
#nBigMatrix = 20, minSize = -5, maxSize = 20, cov = FALSE, envMat = NULL)
#Flist <- makeListFmatrix(fecObjList = fvlist[1:10], nBigMatrix = 20, minSize = -5, 
#maxSize = 20, cov = FALSE, envMat = NULL)
#res <- getIPMoutput(TmatrixList = Tlist, targetSize = 10, FmatrixList = Flist)
#names(res)

# Plot out with different colours for different rows:
#par(mfrow=c(2,2),bty="l",pty="s")
#plot(Tlist[[1]]@meshpoints,res$LE[1, ],xlab="Continuous (e.g. size) stage", 
#ylab="Life expectancy", type="l", 
#	ylim=range(res$LE,na.rm=TRUE))
#for (j in 1:nrow(res$LE)) points(Tlist[[1]]@meshpoints,res$LE[j,],col=j, type="l")

#plot(Tlist[[1]]@meshpoints,res$pTime[1, ], xlab = "Continuous (e.g. size) stage", 
#ylab = "Passage time", type = "l", ylim = range(res$pTime, na.rm=TRUE))
#for (j in 1:nrow(res$pTime)) points(Tlist[[1]]@meshpoints, res$pTime[j, ], col = j, 
#type = "l")

#plot(Tlist[[1]]@meshpoints, Re(res$stableStage[1, ]), xlab = "Continuous (e.g. size) stage", 
#	ylab = "Stable stage distribution", type = "l", 
#	ylim = range(Re(res$stableStage), na.rm=TRUE))
#for (j in 1:nrow(res$stableStage)) points(Tlist[[1]]@meshpoints, 
#	Re(res$stableStage[j, ]),col = j, type = "l")

#hist(res$lambda, xlab = expression(lambda), ylab = "", main = "", col = "grey")
}
