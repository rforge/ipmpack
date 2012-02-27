\name{StochPassageTime}
\alias{StochPassageTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate passage time in a discretely varying environment
}
\description{
  Estimates the time in units of the chosen time-steps (see convertIncrement()) that it will take to
  reach a chosen size for the first time conditional on surviving from
  each of the meshpoints of the IPM for each starting environment. 
}
\usage{
StochPassageTime(chosen.size, IPM.matrix, envMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chosen.size}{
    numeric, target size for which passage time is desired
}
  \item{IPM.matrix}{
    object of class IPM.matrix describing survival related transitions only
}
  \item{envMatrix}{
    object of class envMatrix describing transitions between discrete
    environmental states
}
}
\details{
  Passage time for values exactly equal to the target size are one year, because of how the conditionals are
  framed. Values slightly less than the target size may on average take
  longer due to variance in growth, mortality, leading to
  discontinuities in the pattern of passage time over age. Passage time
  from values > than the chosen.size should be ignored (space to the right of the red vertical line in example below) 
}
\value{
  vector containing passage times across size in the first discrete
  environmental state, then in the second, etc. 
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. 

Metcalf, Horvitz, Tuljapurkar, Clark. 2009. A time to grow and a time to die: a new way to analyze the dynamics of size, light, age and death of
tropical trees. Ecology 90, p2766-2778

For bias in this estimation where variance in growth is small relative to the size range: Zuidema, Jongejans, Chien, During & Schieving. 2010. Integral Projection Models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
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
Tmatrix <-create.compound.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),
growObj=makeGrowthObj(dff, explanatoryVariables="size+size2+covariate",),
survObj=makeSurvObj(dff,explanatoryVariables="size+covariate"), 
envMatrix=makeEnvObj(dff))
targetSize=8

passage=StochPassageTime(targetSize,Tmatrix,makeEnvObj(dff))

plot(passage,ylab="Passage time",xlab="Continuous stage",type="l",col="dark gray",ylim=c(0,max(passage)))
abline(v=targetSize,col="red")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
