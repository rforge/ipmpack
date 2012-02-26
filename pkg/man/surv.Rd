\name{surv}
\alias{surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival
}
\description{
Function to predict the probability of surviving at a given size given a survival object
}
\usage{
surv(size, cov, survObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
    a numeric vector of current sizes
}
  \item{cov}{
    a numeric covariate level (numeric of length 1)
}
  \item{survObj}{
    a survObj 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a vector of length size with values between 0 and 1
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{growth}}, \code{\link{growSurv}}, \code{\link{surv-methods}}
}
\examples{
dff<-generateData()
sv1<-makeSurvObj(dff)
surv(1:50,1,sv1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
