\name{makeListTmatrix}
\alias{makeListTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Builds list of T matrices from list growth and survival objects.
}
\description{
  Uses posteriors (or other) lists of growth and survival objects to create a list of T matrices.
}
\usage{
makeListTmatrix(growObjList, survObjList, nBigMatrix, 
	minSize, maxSize, cov = FALSE, envMat = NULL, integrateType = "midpoint",
		correction = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growObjList}{
    list of growth objects.
  }
  \item{survObjList}{
    list of survival objects.
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
    level of the covariate.
  }
  \item{envMat}{
    environmental matrix - defaults to NULL.
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
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListTmatrix}}, \code{\link{makeListIPMs}}
}
\examples{
dff <- generateData()

svlist <- makePostSurvivalObjs(dataf = dff, explanatoryVariables="size+size2",meanB = rep(0, 3), 
	varB = rep(1e+10), nitt = 5000)
grlist <- makePostGrowthObjs(dataf = dff,explanatoryVariables="size+size2",burnin=100, nitt = 500)

# Only first 10 for speed
Tlist <- makeListTmatrix(growObjList = grlist[1:10], 
	survObjList = svlist[1:10], nBigMatrix = 20, minSize = -5, 
		maxSize = 20, cov = FALSE, envMat = NULL)

# Plots four T matrices at random:
par(mfrow=c(2,2))
image(Tlist[[1]]@meshpoints, Tlist[[1]]@meshpoints, 
	t(Tlist[[1]]), xlab = "Continuous (e.g. size) stage at time t", 
		ylab = "Continuous (e.g. size) stage at time t+1")

image(Tlist[[2]]@meshpoints, Tlist[[2]]@meshpoints, 
	t(Tlist[[2]]), xlab = "Continuous (e.g. size) stage at time t", 
		ylab = "Continuous (e.g. size) stage at time t+1")

image(Tlist[[3]]@meshpoints, Tlist[[3]]@meshpoints, 
	t(Tlist[[3]]), xlab = "Continuous (e.g. size) stage at t", 
		ylab = "Continuous (e.g. size) stage at t+1")

image(Tlist[[4]]@meshpoints, Tlist[[4]]@meshpoints, 
	t(Tlist[[4]]), xlab = "Continuous (e.g. size) stage at time t", 
		ylab = "Continuous (e.g. size) stage at time t+1")
}
