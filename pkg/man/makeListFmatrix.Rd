\name{makeListFmatrix}
\alias{makeListFmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds list of F (sexual reproduction) matrices from lists of growth, survival and fecundity objects.
}
\description{
Uses posteriors (or other) lists of growth, fecundity and survival objects to create a list of F (sexual reproduction) matrices.
}
\usage{
makeListFmatrix(fecObjList, nBigMatrix, minSize, 
	maxSize, cov = FALSE, envMat = NULL, integrateType = "midpoint",
		correction = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{fecObjList}{
	list of fecundity objects.
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
  \item{cov}{
    level of covariate.
  }
  \item{envMat}{
	environmental matrix, defaults to NULL.
}
  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density
    function); other option is "cumul" (which uses the cumulative density function).
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
list of Fmatrices. 
}
\references{
Yang, Jongejans, Yang & Bishop. 2011. The effect of consumers and mutualists of Vaccinum membranaceum at Mount St. Helens: dependence on successional context. PLoS One 10, p1-11.

}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListTmatrix}}, \code{\link{makeListIPMs}}
}
\examples{
dff <- generateData()
fvlist <- makePostFecObjs(dataf = dff, nitt = 5000)

# Only first 10 for speed
Flist <- makeListFmatrix(fecObjList = fvlist[1:10], 
	nBigMatrix = 20, minSize = -5, maxSize = 20, cov = FALSE, 
		envMat = NULL)

# Plots four F matrices at random:
par(mfrow=c(2,2))
image(Flist[[1]]@meshpoints, Flist[[1]]@meshpoints, 
	t(Flist[[1]]), xlab = "Continuous (e.g. size) stage at time t", 
		ylab = "Continuous (e.g. size) stage at time t+1")

image(Flist[[2]]@meshpoints, Flist[[2]]@meshpoints, 
	t(Flist[[2]]), xlab = "Continuous (e.g. size) stage at time t", 
	ylab = "Continuous (e.g. size) stage at time t+1")

image(Flist[[3]]@meshpoints, Flist[[3]]@meshpoints, 
	t(Flist[[3]]), xlab = "Continuous (e.g. size) stage at t", 
	ylab = "Continuous (e.g. size) stage at t+1")

image(Flist[[4]]@meshpoints, Flist[[4]]@meshpoints, 
	t(Flist[[4]]), xlab = "Continuous (e.g. size) stage at time t", 
	ylab = "Continuous (e.g. size) stage at time t+1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
