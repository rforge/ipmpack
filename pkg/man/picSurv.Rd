\name{picSurv}
\alias{picSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Function to make pictures of survival
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
picSurv(dataf, survObj, ncuts = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
%%     ~~Describe \code{dataf} here~~
}
  \item{survObj}{
%%     ~~Describe \code{survObj} here~~
}
  \item{ncuts}{
%%     ~~Describe \code{ncuts} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
dff<-generateData()
sv1 <- makeSurvObjGeneral(dff)
picSurv(dff,sv1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
