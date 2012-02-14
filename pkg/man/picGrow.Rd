\name{picGrow}
\alias{picGrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
picGrow(dataf, growObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
%%     ~~Describe \code{dataf} here~~
}
  \item{growObj}{
%%     ~~Describe \code{growObj} here~~
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
function (dataf, growObj) 
{
    plot(dataf$size, dataf$sizel, pch = 19, xlab = "Size at t", 
        ylab = "Size at t+1", main = "Growth")
    if (length(growObj@fit$coefficients) > 3) {
        dataf$covariate <- as.factor(dataf$covariate)
        levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
        ud <- unique(dataf$covariate)
        ud <- ud[!is.na(ud)]
        for (k in 1:length(ud)) {
            tp <- dataf$covariate == ud[k]
            points(dataf$size[tp], dataf$sizel[tp], pch = 19, 
                col = k)
        }
        ud <- as.factor(ud)
    }
    else {
        ud <- 0
    }
    print(ud)
    sizes <- dataf$size[!is.na(dataf$size)]
    sizes <- sizes[order(sizes)]
    for (k in 1:length(ud)) {
        newd <- data.frame(size = sizes, size2 = sizes^2, logsize = log(size), 
            logsize2 = (log(size))^2, covariate = rep(ud[k], 
                length(sizes)))
        pred.size <- predict(growObj@fit, newd, type = "response")
        if (class(growObj) == "growthObj" | class(growObj) == 
            "growthObj.declinevar") {
            points(sizes, pred.size, type = "l", col = k + 2)
        }
        else {
            if (class(growObj) == "growthObj.incr") {
                points(sizes, sizes + pred.size, type = "l", 
                  col = k + 2)
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
