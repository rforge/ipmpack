\name{SizeToAge}
\alias{SizeToAge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for estimating size-age relationships
}
\description{
  Uses a T matrix and a starting size to estimate time to a range of target
sizes, from which a size to age pattern can be determined. 
}
\usage{
SizeToAge(Tmatrix, starting.size, target.size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrix}{
    object of class IPM.matrix describing growth survival transitions. 
}
  \item{starting.size}{
    numeric, size at age 1, or size from which age-size relationship is desired
}
  \item{target.size}{
    vector of size(s) for which first passage time is required
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{timeInYears}{Time taken to reach target sizes in unit of the time-step
  of the model}
  \item{target.size}{Vector of target sizes}
  \item{starting.size}{Size at age 1}

}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p110-132

Metcalf, Horvitz, Tuljapurkar & Clark. 2009 A time to grow and a time to die: a 
new way to analyze the dynamics of size, light, age and death of
tropical trees. Ecology 90, p2766-2778

Cochran & Ellner. 1992. Simple methods for calculating age-based life history parameters for stage-structured populations. Ecological Monographs 62, p345-364

See for bias in this estimation where variance in growth
is small relative to the size range: Zuidema, Jongejans, Chien, During & Schieving. 2010. Integral Projection Models for trees: a new
parameterization method and a validation of model output. Journal of
Ecology 98, p345-355 

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
Tmatrix <-create.IPM.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
target.size=8
starting.size=0
SizeToAge(Tmatrix, starting.size, target.size)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
