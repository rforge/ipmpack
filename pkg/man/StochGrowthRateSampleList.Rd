\name{StochGrowthRateSampleList}
\alias{StochGrowthRateSampleList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimating the stochastic population growth rate
}
\description{
  Function to estimate the stochastic growth rate lambda_s by iteration; operates by sampling a list of IPMs
}
\usage{
StochGrowthRateSampleList(list.IPM.matrix, n.runin, Tmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.IPM.matrix}{
 a list of IPM.matrix objects corresponding to possible states of the environment
}
  \item{n.runin}{
    numeric, the size of the burnin
}
  \item{Tmax}{
    numeric, the total samples desired
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric converging on high enough log lambda_s for sufficient Tmax; note that if the population size declines to zero, this may return NAs because of logging

}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p452-502

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListIPMs}}
}
\examples{
dff<-generateData()
IPMlist <- makeListIPMs(dataf = dff, nBigMatrix = 10, minSize = -5, maxSize = 25,fec.constants=1000)
StochGrowthRateSampleList(list.IPM.matrix = IPMlist,n.runin = 100, Tmax = 5000) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
