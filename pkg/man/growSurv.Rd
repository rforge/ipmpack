\name{growSurv}
\alias{growSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function combining growth and survival
}
\description{
Function to predict the pdf of size at t+1 given size at t and survival
probability given a growth and survival object
}
\usage{
growSurv(size, sizeNext, cov, growthObj, survObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
    vector of current size(s)
}
  \item{sizeNext}{
    vector of future size(s)
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
\code{create.IPM.Tmatrix} and variants there-of apply \code{outer} to this function to 
efficiently obtain the IPM T matrix. 
}
\value{
numeric defining the pdf
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
Code developed following Mark Rees, Dylan Childs, Karen Rose. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{growth}},  \code{\link{surv}}
}
\examples{
dff<-generateData()
gr1<-makeGrowthObj(dff)
sv1 <- makeSurvObj(dff)
growSurv(1,1:5,1,gr1,sv1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
