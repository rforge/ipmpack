\name{SizeToAge}
\alias{SizeToAge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for estimating size-age relatinoship
}
\description{
  Uses a Tmatrix and a starting size to estimate time to a range of target
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
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PassageTime}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
