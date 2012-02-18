\name{varPassageTime}
\alias{varPassageTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variation in passage time
}
\description{
Function to take a Tmatrix (either compound or not) and estimate variance in passage time to a chosen size
}
\usage{
varPassageTime(chosen.size, IPM.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chosen.size}{
The size of interest
}
  \item{IPM.matrix}{
The Tmatrix (compound or not)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric vector corresponding to variance in passage time from each of the meshpoints in the IPM (so both size, and if a compound matrix, size from different environments). 
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. Sinauer. p119
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PassageTime}}, \code{\link{create.IPM.Tmatrix}} 
}
\examples{
dff <- generateData()
start.sizes <- rnorm(100,10,1)
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff), survObj=makeSurvObj(dff))
vP <- varPassageTime(5,Tmatrix)
plot(Tmatrix@meshpoints,vP, type="l", xlab="Meshpoints", ylab="Variance in passage time")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
