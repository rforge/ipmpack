\name{makeEnvObj}
\alias{makeEnvObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to build environmental transition objects}
\description{
  Function that takes vectors of discrete environmental states at t and
  t+1 and builds a transition matrix from this.
}
\usage{
makeEnvObj(dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `covariate'
    and `covariatel' indicating environmental covariate values at t,
    and at t+1, respectively;
    these must take values of sequential integers, starting at
    `1'}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  an object of class env.matrix
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeFecObj}}, \code{\link{makeSurvObj}}, \code{\link{makeEnvObj}}

}
\examples{
dff <- generateData()
env <- makeEnvObj(dff)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
