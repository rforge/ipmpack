\name{makeGrowthObjIncrTrunc}
\alias{makeGrowthObjIncrTrunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeGrowthObjIncrTrunc(dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
%%     ~~Describe \code{dataf} here~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dataf) 
{
    require(censReg)
    dataf$size2 <- dataf$size^2
    if (length(dataf$incr) == 0) 
        dataf$incr <- dataf$sizenext - dataf$size
    fit <- censReg(incr ~ logsize + logsize2, data = dataf, left = 0)
    gr1 <- new("growthObj.truncincr")
    gr1@fit <- fit$estimate
    return(gr1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
