\name{getIPMoutputDirect}
\alias{getIPMoutputDirect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get IPM output from a list of survival, growth, and (facultatively)
fecundity objects
}
\description{
Get synthetic values including life expectancy, passage time, and if a
fecundity matrix is available, lambda, stable size, etc., for e.g.,
posterior distribution of survival and growth objects. This is essentially identical,
but can be a slightly swifter approach than getIPMoutput since matrices are not stored. It 
also allows compound matrices to be produced following levels of a discrete covariate (e.g., habitat) (using cov=TRUE). 
}
\usage{
getIPMoutputDirect(survObjList, growObjList, target.size = c(), 
	nBigMatrix, minSize, maxSize, cov = FALSE, 
	fecObjList = NULL, env.mat = NULL, n.size.to.age = 0, 
	size.start = 10, integrate.type = "midpoint",correction = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survObjList}{
	list of survival objects
}
  \item{growObjList}{
	list of growth objects
}
  \item{target.size}{
    Target size for passage time. If none is provided defaults to the median of the IPM meshpoints
}
  \item{nBigMatrix}{
	number of meshpoints
}
  \item{minSize}{
minimum size
}
  \item{maxSize}{
maximum size
}
  \item{cov}{
boolean identifying if a discrete covariate should be fitted; if cov=TRUE, 
an environmental matrix MUST be supplied. 
}
  \item{fecObjList}{
list of fecundity objects (facultative)
}
  \item{env.mat}{
    environmental matrix (facultative unless cov=TRUE)
}
  \item{n.size.to.age}{
    numeric describing how many size to age defined (0 - 100s)
}
  \item{size.start}{
    starting size for size to age
}
  \item{integrate.type}{
    integration type, defaults to "midpoint" (which uses prob. density
    function); other option is "cumul" (which uses the cumul. density function)
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
  \item{ptime}{matrix of passage times to the target.size from each of
    the meshpoints (columns) and for each element in the T matrix list (columns)}
  \item{lambda}{vector of lambda corresponding to value obtained
    combining each element of the list of T matrices with the
    corresponding element in the list of F matrices; if no F matrix list is provided, this is a vector of NAs}
  \item{stable.size}{matrix of stable size distributions rows corresponding to values obtained
    combining each element of the list of T matrices with the
    corresponding element in the list of F matrices; if no F matrix list
    is provided, this is a matrix of NAs}
  \item{meshpoints} {matrix meshpoints}
  \item{resSize} {matrix providing target sizes for size to age estimate (assuming age=1
    at size.start), of length n.size.to.age space equally between the smallest and largest meshpoints}
  \item{resAge} {matrix providing time in time-steps to get to resSize,
    rows corresonding sequential elements in the list of T matrices} 
  \item{surv.par} {matrix of posteriors of survival parameters (columns,
    each sample is a row)}
  \item{grow.par} {matrix of posteriors of growth parameters (in columns)} 
}
\references{
For midpoint: Zuidema, Jongejans, Chien, Duiring, Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98 p345-355}
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
dff <- generateData()

#use package MCMCglmm to get Bayesian posteriors - note that sometimes
#!this sometimes later leads to error messages with glm - so you might have
#to re-start R to run other parts of the package! 

svlist <- makePostSurvivalObjs(dataf=dff, meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
grlist <- makePostGrowthObjs(dataf=dff,  nitt=5000)
fvlist <- makePostFecObjs(dataf=dff, nitt=5000)

#Just use first 10 for speed
res <- getIPMoutputDirect(survObjList=svlist[1:10], growObjList=grlist[1:10], target.size = 10,
nBigMatrix=20, minSize=-2, maxSize=15, cov = FALSE, fecObjList =fvlist[1:10],
env.mat = NULL, n.size.to.age = 5, size.start = -1)

names(res)



## plot out with different colours for different rows
par(mfrow=c(2,2),bty="l",pty="s")
plot(res$meshpoints,res$LE[1,],xlab="Size", ylab="Life expectancy", type="l", 
	ylim=range(res$LE,na.rm=TRUE))
for (j in 1:nrow(res$LE)) points(res$meshpoints,res$LE[j,],col=j, type="l")

plot(res$meshpoints,res$ptime[1,],xlab="Size", ylab="Passage Time", type="l", 
	ylim=range(res$ptime,na.rm=TRUE))
for (j in 1:nrow(res$ptime)) points(res$meshpoints,res$ptime[j,],col=j, type="l")

plot(res$meshpoints,Re(res$stable.size[1,]),xlab="Size", ylab="Stable size dist", type="l", 
	ylim=range(Re(res$stable.size),na.rm=TRUE))
for (j in 1:nrow(res$stable.size)) points(res$meshpoints,Re(res$stable.size[j,]),col=j, type="l")

hist(res$lambda,xlab=expression(lambda),ylab="", main="",col="grey")


# same but for compound matrix #
res <- getIPMoutputDirect(survObjList=svlist[1:10], growObjList=grlist[1:10], target.size = 10,
nBigMatrix=20, minSize=-2, maxSize=15, cov = TRUE, fecObjList =fvlist[1:10],
env.mat = makeEnvObj(dff), n.size.to.age = 5, size.start = -1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
