\name{varPassageTime}
\alias{varPassageTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
variation in passage time
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varPassageTime(chosen.size, IPM.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chosen.size}{
%%     ~~Describe \code{chosen.size} here~~
}
  \item{IPM.matrix}{
%%     ~~Describe \code{IPM.matrix} here~~
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
structure(function (chosen.size, IPM.matrix) 
standardGeneric("varPassageTime"), generic = structure("varPassageTime", package = "IPMpack"), package = "IPMpack", group = list(), valueClass = character(0), signature = c("chosen.size", 
"IPM.matrix"), default = `NULL`, skeleton = function (chosen.size, 
    IPM.matrix) 
stop("invalid call in method dispatch to \"varPassageTime\" (no default method)", 
    domain = NA)(chosen.size, IPM.matrix), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
