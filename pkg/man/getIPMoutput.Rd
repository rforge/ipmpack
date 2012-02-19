\name{getIPMoutput}
\alias{getIPMoutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get IPM output from a list of T matrices (only survival and size information) or full IPMs (T matrices and F matrices; the latter include sexual reproduction information)
}
\description{
Get synthetic values including life expectancy, passage time, and, if a
fecundity matrix is available, lambda, stable size, reproductive output, etc. For example, posterior distribution of T matrices, or values across environments, etc. 
}
\usage{
getIPMoutput(Tmatrixlist, target.size = c(), 
	Fmatrixlist = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrixlist}{
    List of survival-size (T) matrices for which summary statistics desired
}
  \item{target.size}{
    Target size for passage time - if none is provided defaults to the median of the IPM meshpoints
}
  \item{Fmatrixlist}{
    List of fecundity (F) matrices corresponding to the list of Tmatrices, facultative. When this information is included, lambda, and stable size distribution will also be provided
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{LE}{matrix of life expectancies, columsn correspond to
  meshpoints, rows corresponding to each element of the list of Tmatrices}
  \item{ptime}{matrix of passage times to the target.size from each of
    the meshpoints (columns) and for each element in the Tmatrix list
    (columns)}
  \item{lambda}{vector of lambda corresonding to value obtained
    combining each element of the list of Tmatrices with the
    corresponding element in the list of Fmatrices; if no Fmatrix list
    is provided, this is a vector of NAs}
  \item{stable.size}{matrix of stable size distributions rows corresonding to values obtained
    combining each element of the list of Tmatrices with the
    corresponding element in the list of Fmatrices; if no Fmatrix list
    is provided, this is a matrix of NAs}
  
}
\references{
Zuidema, Jongejans, Chien, During & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98 p345-355


}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{help}} {getIPMoutputDirect.Rd}
}
\examples{
dff <- generateData()

#use package MCMCglmm to get Bayesian posteriors - note that sometimes
#!this may lead later to error messages with glm. User might have
#to re-start R to run other parts of the package 
svlist <- makePostSurvivalObjs(dataf=dff, explanatoryVariables="size",
meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
grlist <- makePostGrowthObjs(dataf=dff,  nitt=5000)
fvlist <- makePostFecObjs(dataf=dff, nitt=5000)

#use only first 10 of each list, for speed
Tlist <- makeListTmatrix(growObjList=grlist[1:10], survObjList=svlist[1:10],
n.big.matrix=20, minSize=-5, maxSize=20, cov = FALSE, env.mat = NULL)

Flist <- makeListFmatrix(growObjList=grlist[1:10], survObjList=svlist[1:10],fecObjList=fvlist[1:10],
n.big.matrix=20, minSize=-5, maxSize=20, cov = FALSE, env.mat = NULL)

res <- getIPMoutput(Tmatrixlist=Tlist, target.size = 10, Fmatrixlist = Flist)
names(res)


## plot out with different colours for different rows
par(mfrow=c(2,2),bty="l",pty="s")
plot(Tlist[[1]]@meshpoints,res$LE[1,],xlab="Size", ylab="Life expectancy", type="l", 
	ylim=range(res$LE,na.rm=TRUE))
for (j in 1:nrow(res$LE)) points(Tlist[[1]]@meshpoints,res$LE[j,],col=j, type="l")

plot(Tlist[[1]]@meshpoints,res$ptime[1,],xlab="Size", ylab="Passage Time", type="l", 
	ylim=range(res$ptime,na.rm=TRUE))
for (j in 1:nrow(res$ptime)) points(Tlist[[1]]@meshpoints,res$ptime[j,],col=j, type="l")

plot(Tlist[[1]]@meshpoints,res$stable.size[1,],xlab="Size", ylab="Stable size dist", type="l", 
	ylim=range(res$stable.size,na.rm=TRUE))
for (j in 1:nrow(res$stable.size)) points(Tlist[[1]]@meshpoints,res$stable.size[j,],col=j, type="l")

hist(res$lambda,xlab=expression(lambda),ylab="", main="",col="grey")

#TODO - add warning to function

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
