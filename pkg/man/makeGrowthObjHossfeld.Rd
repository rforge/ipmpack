\name{makeGrowthObjHossfeld}
\alias{makeGrowthObjHossfeld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeGrowthObjHossfeld(dataf)
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
function (dataf) 
{
    if (length(dataf$incr) == 0) 
        dataf$incr <- dataf$sizel - dataf$size
    dataf$incr[dataf$incr < 0] <- 0
    tmp <- optim(c(1, 1, 1), wrapHossfeld, dataf = dataf, method = "Nelder-Mead")
    print(tmp$convergence)
    gr1 <- new("growthObj.Hossfeld")
    gr1@paras <- tmp$par
    resids <- Hossfeld(dataf$size, tmp$par) - dataf$incr
    gr1@sd <- sd(resids, na.rm = T)
    return(gr1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
