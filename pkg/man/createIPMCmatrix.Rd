\name{createIPMCmatrix}
\alias{createIPMCmatrix}
%- Also NEED an '\alias' for EACH other topic documented here. 
\title{
Builds C matrices.
}
\description{
Uses clonality objects to construct a matrix defining per-capita
contribution to clonal stages (e.g., propagules [seed, spore],
seedlings, calves) by clonal reproduction. Currently only pre-census
clonal reproduction can be handled. 
}
\usage{
createIPMCmatrix(clonalObj, nEnvClass = 1, nBigMatrix = 50, minSize = -1, 
	maxSize = 50, chosenCov = data.frame(covariate=1), integrateType="midpoint", correction="none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clonalObj}{
    clonal reproduction object; currently essentially identical to a
    fecundity reproduction object
  }
  \item{nEnvClass}{
    numeric, number of environmental classes, always = 1 for
    non-compound matrices.
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
  \item{chosenCov}{
   data-frame indicating level of the discrete covariate, or 
    range of values where multiple covariates are modeled.
  }
  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density
    function); other option is "cumul" (which uses the cumulative density function)
  }
  \item{correction}{
  	correction type, defaults to \code{none}. The first option is \code{constant} which will multiply every column of the IPM by a constant sufficient to adjust values to those predicted for total fertility at that size. The second option is \code{discretizeExtremes} which will place all transitions to sizes smaller than \code{minSize} into the smallest bin, and transitions to sizes larger than \code{maxSize} into the largest bin.  
   }
 
}
\value{
  an object of class IPMmatrix of dimensions nBigMatrix or nBigMatrix+nDiscrete classes (defined by clonalObj@offspringSplitter-1).
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
\code{\link{createIPMPmatrix}},\code{\link{createIPMFmatrix}}

}
\examples{
# Data with only continuous stage and one habitat
dff <- generateData()
dff$fec[dff$fec==0] <- NA
cv1 <- makeClonalObj(dff, Formula = log(fec)~size, Transform = "log")
Cmatrix <- createIPMCmatrix(clonalObj = cv1, nBigMatrix = 20, 
	minSize = min(dff$size, na.rm = TRUE), maxSize = max(dff$size, na.rm = TRUE))

slotNames(Cmatrix)

image(Cmatrix@meshpoints, Cmatrix@meshpoints, t(Cmatrix), 
	xlab = "Continuous (e.g. size) stage at t", 
		ylab = "Continous (e.g. size) stage at t+1")
}
