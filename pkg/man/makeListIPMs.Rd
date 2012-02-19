\name{makeListIPMs}
\alias{makeListIPMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make list IPMs where there is a discretely covarying covariate
}
\description{
Wrapper function to build the IPM corresponding to every level of the
discrete covariate, and return a list of these. 
}
\usage{
makeListIPMs(dataf, nBigMatrix = 10, minSize = -2, maxSize = 10,
	integrate.type="midpoint", correction="none",
	explSurv="size+size2+covariate",explGrow="size+size2+covariate", 
	regType="constantVar",responseType="sizeNext",explFec="size",fec.constants=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size', `sizeNext', 'surv', 'fec',
    'covariate', 'covariatel'; 
    and 'age' indicating which individuals
    are seedlings for identifying the mean and variance in seedling
      size. }
 \item{nBigMatrix}{
	number of bins in size
	}
  \item{minSize}{
	smallest size
	}
  \item{maxSize}{
	largest size
	}
  \item{integrate.type}{
	integration type
	}
  \item{correction}{
	correction
	}
 \item{explSurv}{
	explanatory variables used in the survival model
	}
 \item{explGrow}{
	explanatory variables used in the growth model
	}
 \item{regType}{
	regression Type for growth (normal density function, truncated....)
	}
 \item{responseType}{
	regression response variable for growth (sizeNext, size...)
	}
 \item{explFec}{
	explanatory variables used in the fecundity
	}
 \item{fec.constants}{
	constant multipliers for the fecundity model
	}
}
\details{
Currently uses makeSurvObj with a discrete covariate 
}
\value{
  list of matrices corresponding to covariates, in order
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
dff <- generateData()
tmp <- makeListIPMs(dff)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
