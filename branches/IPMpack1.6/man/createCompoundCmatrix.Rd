\name{createCompoundCmatrix}
\alias{createCompoundCmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds a compound C matrix.
}
\description{
  Uses clonality object, and environmental transition objects to
  construct a matrix defining probabilities for transitions between
  sizes due to clonal reporudction given both a continuous state
  (e.g. size) and environmental state, as well as a discrete stage if
  necessary. Currently only pre-census clonal reproduction can be handled.  
}
\usage{
createCompoundCmatrix(nEnvClass = 2, nBigMatrix = 50,
minSize = -1, maxSize = 50, envMatrix, clonalObj, integrateType, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nEnvClass}{
    numeric, number of environmental classes, defaults to 2.
  }
  \item{nBigMatrix}{
    numeric, number of size bins in the P matrix, defaults to 50.
  }
  \item{minSize}{
    numeric, minimum size of the P matrix, defaults to -1.
  }
  \item{maxSize}{
    numeric, maximum size of the P matrix, defaults to 50.
  }
  \item{envMatrix}{
    envMatrix object defining transitions between environmental states
    for each size.
  }
  \item{clonalObj}{
    clonality object.
  }
  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density
    function); other option is "cumul" (which uses the cumulative density function).
  }
  \item{correction}{
 correction type, defaults to \code{none}. The first option is \code{constant} which will multiply every column of the IPM by a constant sufficient to adjust values to those predicted for total fertility at that size. The second option is \code{discretizeExtremes} which will place all transitions to sizes smaller than \code{minSize} into the smallest bin, and transitions to sizes larger than \code{maxSize} into the largest bin.    }
}
\value{
  an object of class IPMmatrix with dimensions
  nBigMatrix*nEnvClass, or if discrete transitions exist
  (nBigMatrix+nDisc)*nEnvClass

}
\references{
For information on C matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.

For midpoint: Zuidema, Jongejans, Chien, During & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355.

For multiple-vital rate integration on fecundity: Yang, Jongejans, Yang & Bishop. 2011. The effect of consumers and mutualists of Vaccinum membranaceum at Mount St. Helens: dependence on successional context. PLoS One 10, p1-11.

}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createCompoundPmatrix}},\code{\link{createIPMCmatrix}}
}
\examples{
## See createCompoundFmatrix for examples
}
