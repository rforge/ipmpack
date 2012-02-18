\name{timeToSize}
\alias{timeToSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to project how long it takes to get from a starting distribution to an final size
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
timeToSize(starting.sizes, IPM, end.size, starting.env = 1, maxt = 100, prop.reach = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{starting.sizes}{
vector of starting sizes (in any order)
}
  \item{IPM}{
the IPM (just a Tmatrix)
}
  \item{end.size}{
the end size
}
  \item{starting.env}{
vector of starting env, same length as starting.sizes, or length=1
}
  \item{maxt}{
the max number of time-steps tested
}
  \item{prop.reach}{
the proportion of the starting pop that have to be > than the end.size for it to count
}
}
\details{
Plots and returned values of survivorship from preliminary runs will give a notion of how low this has to be
}
\value{
 \item{ts.dist}{ the time-series of size distribution}
 \item{time.reach }{the time for n.reach to be at sizes > end.size }
 \item{survivorship}{survivorship over the course of the time elapsed for that pop}
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. Sinauer. 
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
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
start.sizes <- rnorm(100,2,1)
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObjGeneral(dff), survObj=makeSurvObjGeneral(dff))

#CURRENTLY NOT WORKING VERY WELL
rc <- timeToSize(starting.sizes=start.sizes, IPM=Tmatrix, end.size=10, 
starting.env = 1, maxt = 1000, prop.reach = 0.001)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
