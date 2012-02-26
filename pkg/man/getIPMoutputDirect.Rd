\name{getIPMoutputDirect}
\alias{getIPMoutputDirect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gets IPM output from a list of survival, growth, and (facultatively)
fecundity objects.
}
\description{
Get synthetic values including life expectancy, passage time, and if a
fecundity matrix is available, population growth rate (lambda), stable stage disctribution, etc. This is essentially identical to getIPMoutput(),but can be a slightly swifter approach since matrices are not stored. It 
also allows compound matrices to be produced following levels of a discrete covariate (e.g., habitat) (using cov=TRUE). 
}
\usage{
getIPMoutputDirect(survObjList, growObjList, targetSize = c(), 
	nBigMatrix, minSize, maxSize, discreteTrans = 1,cov = FALSE, 
	fecObjList = NULL, envMat = NULL, nSizeToAge = 0, 
	sizeStart = 10, integrateType = "midpoint",correction = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survObjList}{
	list of survival objects.
}
  \item{growObjList}{
	list of growth objects.
}
  \item{targetSize}{
    Target size for passage time. If none is provided defaults to the median of the IPM meshpoints.
}
  \item{nBigMatrix}{
	number of meshpoints.
}
  \item{minSize}{
minimum size.
}
  \item{maxSize}{
maximum size.
}
  \item{discreteTrans}{
    object of class discreteTrans, or numeric.
  }
  \item{cov}{
boolean identifying if a discrete covariate should be fitted; if cov=TRUE, 
an environmental matrix MUST be supplied. 
}
  \item{fecObjList}{
list of fecundity objects (facultative).
}
  \item{envMat}{
    environmental matrix (facultative unless cov=TRUE).
}
  \item{nSizeToAge}{
    numeric describing how many size to age defined (0 - 100s).
}
  \item{sizeStart}{
    starting continuous stage (e.g. size) for stage to age.
}
  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density
    function); other option is "cumul" (which uses the cumulative density function)
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will
    multiply every column of the IPM by a constant sufficient to adjust
    values to those predicted for survival at that size. 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{LE}{matrix of life expectancies, columns correspond to
  meshpoints, rows corresponding to each element of the list of T matrices}
  \item{pTime}{matrix of passage times to the targetSize from each of
    the meshpoints (columns) and for each element in the T matrix list (columns).}
  \item{lambda}{vector of population growth rates (lambda) corresponding to value obtained
    combining each element of the list of T matrices with the
    corresponding element in the list of F matrices; if no F matrix list is provided, this is a vector of NAs.}
  \item{stableStage}{matrix of stable size distributions rows corresponding to values obtained
    combining each element of the list of T matrices with the
    corresponding element in the list of F matrices; if no F matrix list
    is provided, this is a matrix of NAs.}
  \item{meshpoints} {matrix meshpoints.}
  \item{resSize} {matrix providing target sizes for size to age estimate (assuming age=1
    at sizeStart), of length nSizeToAge space equally between the smallest and largest meshpoints.}
  \item{resAge} {matrix providing time in time-steps to get to resSize,
    rows corresonding sequential elements in the list of T matrices.} 
  \item{surv.par} {matrix of posteriors of survival parameters (columns,
    each sample is a row).}
  \item{grow.par} {matrix of posteriors of growth parameters (in columns).} 
}
\references{
For midpoint: Zuidema, Jongejans, Chien, During, Schieving. 2010. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Only for continuous stage:
dff <- generateData()

svlist <- makePostSurvivalObjs(dataf = dff, explanatoryVariables = "size+size2", meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
grlist <- makePostGrowthObjs(dataf = dff, explanatoryVariables = "size+size2", meanB = rep(0, 3), varB = rep(1e+10),  nitt=5000)
fvlist <- makePostFecObjs(dataf=dff, explanatoryVariables="size+size2", Transform = "log", nitt=5000)

#Just use first 10 for speed
res <- getIPMoutputDirect(survObjList=svlist[1:10], growObjList=grlist[1:10], targetSize = 10,
nBigMatrix=20, minSize = -2, maxSize = 15, cov = FALSE, fecObjList = fvlist[1:10],
envMat = NULL, nSizeToAge = 5, sizeStart = -1)

names(res)

# Plot out with different colours for different rows:
par(mfrow = c(2,2), bty = "l", pty = "s")
plot(res$meshpoints, res$LE[1, ], xlab = "Continuous (e.g. size) stage", ylab = "Life expectancy", type= "l", ylim = range(res$LE, na.rm=TRUE))
for (j in 1:nrow(res$LE)) points(res$meshpoints,res$LE[j, ],col=j, type = "l")

plot(res$meshpoints, res$pTime[1, ], xlab="Continuous (e.g. Size) stage", ylab = "Passage time", type = "l", ylim = range(res$pTime, na.rm=TRUE))
for (j in 1:nrow(res$pTime)) points(res$meshpoints,res$pTime[j, ],col=j,  type = "l")

plot(res$meshpoints,Re(res$stableStage[1, ]), xlab = "Size", ylab = "Stable stage distribution", type = "l", ylim = range(Re(res$stableStage), na.rm = TRUE))
for (j in 1:nrow(res$stableStage)) points(res$meshpoints, Re(res$stableStage[j, ]), col = j, type = "l")

hist(res$lambda, xlab = expression(lambda), ylab = "", main = "", col = "grey")

## For continuous (e.g. size) and discrete (e.g. seedbank) stages:

svlist <- makePostSurvivalObjs(dataf = dff, explanatoryVariables = "size+size2+covariate", meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
env1 <- makeEnvObj(dff)
res <- getIPMoutputDirect(survObjList=svlist[1:10], growObjList = grlist[1:10], targetSize = 10, nBigMatrix = 20, minSize = -2, maxSize = 15, cov = TRUE, fecObjList =fvlist[1:10], envMat = env1, nSizeToAge = 5, sizeStart = -1)



# same but for matrix with discrete classes #
dff <- generateDataDiscrete()

svlist <- makePostSurvivalObjs(dataf=dff, explanatoryVariables="size+size2", meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
grlist <- makePostGrowthObjs(dataf=dff, explanatoryVariables="size+size2", meanB = rep(0, 3), varB = rep(1e+10),  nitt=5000)
fvlist <- makePostFecObjs(dataf=dff, explanatoryVariables="size+size2", Transform="log",nitt=5000)

res <- getIPMoutputDirect(survObjList=svlist[1:10], growObjList=grlist[1:10], targetSize = 10,
nBigMatrix=20, minSize=-2, maxSize=15, cov = FALSE, fecObjList =fvlist[1:10],
nSizeToAge = 5, sizeStart = -1)





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
