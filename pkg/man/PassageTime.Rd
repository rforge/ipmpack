\name{PassageTime}
\alias{PassageTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to define passage time to a chosen size
}
\description{
  Estimates the time in units of the chosen time-steps that it will take to
  reach a chosen size for the first time conditional on surviving from each of the meshpoints of the IPM
}
\usage{
PassageTime(chosen.size, IPM.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chosen.size}{
    numeric, the target size
}
  \item{IPM.matrix}{
    an IPM.matrix object describing growth-survival transitions (a Tmatrix)
}
}
\details{
  Passage time for values exactly equal to the chosen.size are one year, because of how the conditionals are
  framed. Values slightly less than the target size may on average take
  longer due to variance in growth, mortality, leading to
  discontinuities in the pattern of passage time over age. Passage time
  from values > than the chosen.size should be ignored. 
}
\value{
  a vector of times in the units of the chosen time-steps corresponding to each of the IPM meshpoints
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. Sinauer. 

Metcalf, C.J.E., Horvitz, C.C., Tuljapurkar, S. & Clark, D.A. 2009 A time to grow and a time to die: a 
new way to analyze the dynamics of size, light, age and death of
tropical trees. Ecology 90: 2766-2778

See Zuidema et al. 2010. Integral Projection Models for trees: a new
parameterization method and a validation of model output. Journal of
Ecology 98: 345-355 for bias in this estimation where variance in growth
is small relative to the size range. 
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MeanLifeExpect}}, ~~~
}
\examples{
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size),
maxsize=max(dff$size),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
PassageTime(5,Tmatrix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
