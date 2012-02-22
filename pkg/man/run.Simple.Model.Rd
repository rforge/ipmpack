\name{run.Simple.Model}
\alias{run.Simple.Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to run a simple model using default survival and growth functions. 
}
\description{
This function takes a dataframe, a chosen size, a minimum and maximum size onto which integrate, and a number of bins, and creates a T matrix and plots passage time and life expectancy (does not include fecundity)
}
\usage{
run.Simple.Model(dataf, chosen.size, minSize = c(), maxSize = c(), 
nBigMatrix = 500, do.plot = TRUE, is.log = TRUE,  
integrate.type = "midpoint", correction="none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
  a data-frame with columns size, sizeNext (which is size in the next time
  step), and surv, at a minimum%%     
  }
  \item{chosen.size}{
	the target size for which passage time estimates are desired
}
  \item{minSize}{
    minimum size for the IPM, defaults to a vector of length zero, in which case, 
    it is estimated as either 0.5*min(dataf$size,na.rm=TRUE) if the 
    minimum size is positive, or 1.5*min(dataf$size,na.rm=TRUE) if the minimum size is negative 
    (e.g., reflecting a log scale)
    }
  \item{maxSize}{
    minimum size for the IPM, defaults to a vector of length zero, in which case, 
    it is estimated as 1.5*max(dataf$size,na.rm=TRUE) 
}
  \item{nBigMatrix}{
number of bins used in the IPM
}
  \item{do.plot}{
boolean describing whether plots should be produced
}
  \item{is.log}{
booean describing whether size is on a log scale or not (for plotting)
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
}
\value{
  \item{ptime}{
	passage time to chosen.size from every size in the chosen discrete time step interval; 
	estimation for sizes larger than the chosen.size should be discarded. 
	
  }
  \item{LE}{
	life expectancy  in the chosen discrete time step interval
  }
  \item{Tmatrix}{
	The basic Tmatrix
  }
  \item{growObj}{
	The fitted growth object
  }
  \item{survObj}{
	The fitted survival object
  }

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff<-generateData()
a1<-run.Simple.Model(dataf=dff, chosen.size=3,do.plot=TRUE)
names(a1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
