\name{coerceMatrixIPM}
\alias{coerceMatrixIPM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coerce matrix population model to IPM.matrix bject
}
\description{
Forces MPM (or survival-growth transition matrix, or
similar) into an IPM.matrix object for use with passage time,
life expectancy, etc. functions. 
}
\usage{
coerceMatrixIPM(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{
    A matrix population model
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  an object of class IPM.matrix
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#build simple matrix mode
amat <- matrix(0,3,3)
amat[cbind(2:3,1:2)] <- 0.5
amat[c(3,3)] <- 0.8
amat[1,3] <-5

IPM.amat <- coerceMatrixIPM(amat)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
