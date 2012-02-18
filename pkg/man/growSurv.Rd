\name{growSurv}
\alias{growSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functioning combining growth and survival
}
\description{
Function to predict the pdf of size at t+1 given size at t and survival
probability given a growth and survival object
}
\usage{
growSurv(size, sizenext, cov, growthObj, survObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
    vector of current sizes
}
  \item{sizenext}{
    vector of future sizes
}
  \item{cov}{
    covariate level (numeric of length 1)
}
  \item{growthObj}{
    a growth object
}
\item{survObj}{
  a survival object
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
growth, surv, fertility
}
\examples{
dff<-generateData()
gr1<-makeGrowthObjGeneral(dff)
sv1 <- makeSurvObjGeneral(dff)
growSurv(1,1:5,1,gr1,sv1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
