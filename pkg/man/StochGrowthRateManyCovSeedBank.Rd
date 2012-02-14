\name{StochGrowthRateManyCovSeedBank}
\alias{StochGrowthRateManyCovSeedBank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for estimating lambda_s or invasion rate with many varying covariates
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
StochGrowthRateManyCovSeedBank(covariate, n.runin, Tmax, growthObj, survObj, fecObj, n.big.matrix, minsize, maxsize, n.microsites)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covariate}{
%%     ~~Describe \code{covariate} here~~
}
  \item{n.runin}{
%%     ~~Describe \code{n.runin} here~~
}
  \item{Tmax}{
%%     ~~Describe \code{Tmax} here~~
}
  \item{growthObj}{
%%     ~~Describe \code{growthObj} here~~
}
  \item{survObj}{
%%     ~~Describe \code{survObj} here~~
}
  \item{fecObj}{
%%     ~~Describe \code{fecObj} here~~
}
  \item{n.big.matrix}{
%%     ~~Describe \code{n.big.matrix} here~~
}
  \item{minsize}{
%%     ~~Describe \code{minsize} here~~
}
  \item{maxsize}{
%%     ~~Describe \code{maxsize} here~~
}
  \item{n.microsites}{
%%     ~~Describe \code{n.microsites} here~~
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
function (covariate, n.runin, Tmax, growthObj, survObj, fecObj, 
    n.big.matrix, minsize, maxsize, n.microsites) 
{
    require(MASS)
    nt <- rep(1, n.big.matrix + 1)
    Rt <- rep(NA, Tmax)
    tmp.fecObj <- fecObj
    if (sum(n.microsites) > 0) {
        dd <- TRUE
    }
    else {
        dd <- FALSE
        tmp.fecObj@p.est <- 1
    }
    for (t in 1:Tmax) {
        if (dd) 
            tmp.fecObj@p.est <- min(n.microsites[min(t, length(n.microsites))]/nt[1], 
                1)
        tmp.fecObj@seed.germination <- pgerminate(matrix(covariate, 
            1, length(covariate)), params = tmp.fecObj@germ.params)
        IPM.here <- create.seed.matrix(n.env.class = 1, n.big.matrix = n.big.matrix, 
            minsize = minsize, maxsize = maxsize, chosen.cov = covariate[t, 
                ], growObj = growthObj, survObj = survObj, fecObj = tmp.fecObj)
        nt1 <- IPM.here \%*\% nt
        sum.nt1 <- sum(nt1)
        Rt[t] <- log(sum.nt1)
        nt <- nt1/sum.nt1
    }
    res <- mean(Rt[n.runin:Tmax], na.rm = TRUE)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
