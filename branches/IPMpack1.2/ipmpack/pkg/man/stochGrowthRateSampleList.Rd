\name{stochGrowthRateSampleList}
\alias{stochGrowthRateSampleList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimating the stochastic population growth rate (lambda_s).
}
\description{
  Estimates the stochastic growth rate (lambda_s) by iteration; operates by sampling a list of IPMs.
}
\usage{
stochGrowthRateSampleList(listIPMmatrix, nRunIn, tMax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listIPMmatrix}{
 a list of IPMmatrix objects corresponding to possible states of the environment.
}
  \item{nRunIn}{
    numeric, the size of the burnin,
}
  \item{tMax}{
    numeric, the total samples desired.
}
}
\value{
a numeric converging on high enough log lambda_s for sufficient tMax; note that if the population size declines to zero, this may return NAs because of logging.

}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p452-502.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListIPMs}}
}
\examples{
dff <- generateData()
IPMlist <- makeListIPMs(dataf = dff, nBigMatrix = 10, minSize = -5, 
maxSize = 15,fecConstants=data.frame(1e6), correction="constant")
stochGrowthRateSampleList(listIPMmatrix = IPMlist,nRunIn = 100, tMax = 5000) 
}
