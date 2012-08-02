\name{stochGrowthRateSampleListDD}
\alias{stochGrowthRateSampleListDD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimating the stochastic invasion rate.
}
\description{
 Estimates the stochastic invasion rate in the presence of density dependence in seedling establishment by iteration; operates by sampling a list of IPMs and recalculating the probability of seed establishment at each time-step. }
\usage{
stochGrowthRateSampleListDD(listTmatrix, listFmatrix, nRunIn, tMax,seedList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listTmatrix}{
 a list of IPM T matrices corresponding to possible states of the environment.
}
  \item{listFmatrix}{
 a list of IPM F matrices  corresponding to possible states of the environment.
}
  \item{nRunIn}{
    numeric, the size of the burnin,
}
  \item{tMax}{
    numeric, the total samples desired.
}
  \item{seedList}{
    numeric, a vector of the number of successful recruit corresponding to possible states of the environment. 
}



}
\value{
a numeric converging on high enough log vs for sufficient tMax.

}
\references{
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListFmatrix}}; \code{\link{makeListTmatrix}}; \code{\link{stochGrowthRateSampleList}}

}
\examples{

}
