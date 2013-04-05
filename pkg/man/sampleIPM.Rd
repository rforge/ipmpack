\name{sampleIPM }
\alias{.makeListPmatrix }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Builds list of IPMs or P matrices from list growth, survival, fecundity and discreteTrans objects.
}
\description{
  Uses lists of vital rate objects to create a list of IPM or P matrices.
}
\usage{
.makeListPmatrix (growObjList, survObjList, nBigMatrix, 
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
\code{\link{.makeListPmatrix }}, \code{\link{sampleSequentialIPMs}}
}
\examples{

#TODO: Cory example
}
