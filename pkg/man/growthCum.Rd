\name{growthCum}
\alias{growthCum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model of growth allowing cumulative bin estimation}
\description{
Generic function to predict the cdf of size at t+1 given size at t and a growth object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
growthCum(size, sizeNext, cov, growthObj)
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 a vector of length sizeNext giving the cdf of each value of sizeNext
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
growth
}
\examples{
dff <- generateData()
gr1 <- makeGrowthObjGeneral(dff)
growthCum(1,1:5,1,gr1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
