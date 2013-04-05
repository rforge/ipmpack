\name{.makeListFmatrix }
\alias{.makeListFmatrix }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds list of F (sexual reproduction) matrices from lists of growth, survival and fecundity objects.
}
\description{
Uses posteriors (or other) lists of growth, fecundity and survival objects to create a list of F (sexual reproduction) matrices.
}
\usage{
.makeListFmatrix (fecObjList, nBigMatrix, minSize, 
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
\value{
list of Fmatrices. 
}
\references{
Yang, Jongejans, Yang & Bishop. 2011. The effect of consumers and mutualists of Vaccinum membranaceum at Mount St. Helens: dependence on successional context. PLoS One 10, p1-11.

}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{.makeListPmatrix }}, \code{\link{sampleSequentialIPMs}}
}
\examples{
dff <- generateData()

## TODO: Cory example


}
