\name{makeListIPMs}
\alias{makeListIPMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MakeS list IPMs where there is a discrete covariate.
}
\description{
Wrapper function to build the IPM corresponding to every level of the discrete covariate, and return a list of these. 
}
\usage{
makeListIPMs(dataf, nBigMatrix = 10, minSize = -2, maxSize = 10,	
	integrateType = "midpoint", correction = "none", 
		explSurv = "size+size2+covariate", 
			explGrow = "size+size2+covariate", 
				regType = "constantVar", 
					responseType = "sizeNext", 
						explFec = "size",  Family="gaussian", 
							Transform="none", fecConstants = data.frame(NA))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size', `sizeNext', 'surv', 'fec',
    'covariate', 'covariatel'; and 'age' indicating which individuals are seedlings for identifying the mean and variance in seedling size.
    }
 \item{nBigMatrix}{
	number of bins in size.
	}
  \item{minSize}{
	minimum size.
	}
  \item{maxSize}{
	maximum size.
	}
  \item{integrateType}{
	integration type.
	}
  \item{correction}{
	correction.
	}
 \item{explSurv}{
	explanatory variables used in the survival model.
	}
 \item{explGrow}{
	explanatory variables used in the growth model.
	}
 \item{regType}{
	regression Type for growth (normal density function, truncated, etc).
	}
 \item{responseType}{
	regression response variable for growth (sizeNext, size, etc).
	}
 \item{explFec}{
	explanatory variables used in the fecundity.
	}
  \item{Family}{a character vector containing the names of the families to be used for the glms, e.g., binomial, poisson, etc. Again, these must appear  \emph{in the
      order defined by the alphabetical order of column names 
      containing "fec" or in the order of column names specified in fecNames} }
  \item{Transform}{a character vector containing the names of the
    transforms to be used for the response variables, e.g., log, sqrt, -1, etc. Again, these must appear  \emph{in the order defined by the alphabetical order of column names
      containing "fec" or in the order of column names specified in fecNames} }
 \item{fecConstants}{
	data.frame of constant multipliers for the fecundity model.
	}
}
\details{
Currently uses makeSurvObj with a discrete covariate .
}
\value{
  list of matrices corresponding to covariates, in order.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
dff <- generateData()
IPMlist <- makeListIPMs(dff, Transform="log")
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
