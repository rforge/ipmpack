\name{runSimpleModel}
\alias{runSimpleModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs a simple model using default survival and growth functions. 
}
\description{
Takes a dataframe, a chosen size, a minimum and maximum size onto which integrate, and a number of bins, and creates a P matrix and plots passage time and life expectancy (does not include fecundity).
}
\usage{
runSimpleModel(dataf, chosenSize, minSize = c(), maxSize = c(), 
nBigMatrix = 500, do.plot = TRUE, is.log = TRUE,  
integrateType = "midpoint", correction="none")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
  a dataframe with columns size, sizeNext (which is continuous stage in the next time step), and surv, at a minimum.%%     
  }
  \item{chosenSize}{
	the target size for which passage time estimates are desired.
}
  \item{minSize}{
    minimum size for the IPM, defaults to a vector of length zero, in which case, it is estimated as either 0.5*min(dataf$size,na.rm = TRUE) if the minimum size is positive, or 1.5*min(dataf$size, na.rm = TRUE) if the minimum size is negative (e.g., reflecting a log scale).
    }
  \item{maxSize}{
    minimum size for the IPM, defaults to a vector of length zero, in which case, it is estimated as 1.5*max(dataf$size, na.rm = TRUE). 
}
  \item{nBigMatrix}{
number of bins used in the IPM.
}
  \item{do.plot}{
boolean describing whether plots should be produced.
}
  \item{is.log}{
boolean describing whether continuous stage is on a log scale or not (for plotting).
}

  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density function); other option is "cumul" (which uses the cumulative density function).
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will multiply every column of the IPM by a constant sufficient to adjust values to those predicted for survival at that stage. 
  }
}
\value{
  \item{pTime}{
	passage time to chosenSize from every size in the chosen discrete time step interval; estimation for sizes larger than the chosenSize should be discarded. 
  }
  \item{LE}{
	life expectancy  in the chosen discrete time step interval.
  }
  \item{Pmatrix}{
	The basic P matrix.
  }
  \item{growObj}{
	The fitted growth object.
  }
  \item{survObj}{
	The fitted survival object.
  }

}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
#dff <- generateData()
#a1<-runSimpleModel(dataf=dff, chosenSize=3,do.plot=TRUE)
#names(a1)
}
