\name{growthCum}
\alias{growthCum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{~~ title for growthCum ~~}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
growthCum(size, sizenext, cov, growthObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
%%     ~~Describe \code{size} here~~
}
  \item{sizenext}{
%%     ~~Describe \code{sizenext} here~~
}
  \item{cov}{
%%     ~~Describe \code{cov} here~~
}
  \item{growthObj}{
%%     ~~Describe \code{growthObj} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (size, sizenext, cov, growthObj) 
standardGeneric("growthCum"), generic = structure("growthCum", package = "IPMpack"), package = "IPMpack", group = list(), valueClass = character(0), signature = c("size", 
"sizenext", "cov", "growthObj"), default = `NULL`, skeleton = function (size, 
    sizenext, cov, growthObj) 
stop("invalid call in method dispatch to \"growthCum\" (no default method)", 
    domain = NA)(size, sizenext, cov, growthObj), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
