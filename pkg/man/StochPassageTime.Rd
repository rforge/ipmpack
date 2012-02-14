\name{StochPassageTime}
\alias{StochPassageTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate passage time in a discretely varying environment
}
\description{
  Estimates the time in units of the chosen time-steps that it will take to
  reach a chosen size for the first time conditional on surviving from
  each of the meshpoints of the IPM for each starting environment. 
}
\usage{
StochPassageTime(chosen.size, IPM.matrix, env.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chosen.size}{
    numeric, target size for which passage time is desired
}
  \item{IPM.matrix}{
    object of class IPM.matrix describing survival related transitions only
}
  \item{env.matrix}{
    object of class env.matrix describing transitions between discrete
    environmental states
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
  vector containing passage times across size in the first discrete
  environmental state, then in the second, etc. 
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
\code{\link{PassageTime}}
}
\examples{
dff <- generateData()
Tmatrix <-create.compound.Tmatrix(minsize=min(dff$size),
maxsize=max(dff$size),growObj=makeGrowthObjDiscreteCovariate(dff),
survObj=makeSurvObjDiscreteCovariate(dff), envMatrix=makeEnvObj(dff))
StochPassageTime(10,Tmatrix, makeEnvObj(dff))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
