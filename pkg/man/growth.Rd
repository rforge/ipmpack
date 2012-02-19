\name{growth}
\alias{growth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Growth 
}
\description{
Generic function to predict the pdf (probability density function) of size at t+1 given size at t and a growth object
}
\usage{
growth(size, sizeNext, cov, growthObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
    vector of current sizes
}
  \item{sizeNext}{
    vector of sizes in the next time-step
}
  \item{cov}{
    covariate level in this time-step (numeric of length 1)
}
  \item{growthObj}{
    a growth object
}
}
\details{
  Models defining size at t+1, or growth increment are defined; with various
  underlying statistical models
}
\value{
  a vector of length sizeNext giving the pdf (probability density function) of each value of sizeNext
}
\references{
Easterling, Ellner & Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81 p694-708

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{surv}}, \code{\link{growSurv}}
}
\examples{s
dff <- generateData()
gr1 <- makeGrowthObj(dff)
growth(1,1:5,1,gr1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
