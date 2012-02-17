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
makeListIPMs(dataf, n.big.matrix = 10, minsize = -2, maxsize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size', `sizel', 'surv', 'fec',
    'covariate', 'covariatel'; 
    and 'age' indicating which individuals
    are seedlings for identifying the mean and variance in seedling
      size. }
 \item{n.big.matrix}{
	number of bins in size
	}
  \item{minsize}{
	smallest size
	}
  \item{maxsize}{
	largest size
	}
}
\details{
Currently uses makeSurvObjGeneral with a discrete covariate 
}
\value{
  list of matrices corresponding to covariates, in order
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
