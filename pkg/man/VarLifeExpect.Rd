\name{VarLifeExpect}
\alias{VarLifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating variance in life-expectancy
}
\description{
  Provided with a Tmatrix, defining survival transitions across size,
  provides a vector defining variance in life-expectancy in units of the time-step used, for each of the size bins. }
\usage{
VarLifeExpect(IPM.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPM.matrix}{
    an IPM.matrix object defining survival transitions
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a vector of variance in life-expectancies each corresponding to Tmatrix@meshpoints

}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. Sinauer. p119

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MeanLifeExpect}}, \code{\link{create.IPM.Tmatrix}} 

}
\examples{
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObjGeneral(dff),
survObj=makeSurvObjGeneral(dff))
VarLifeExpect(Tmatrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
