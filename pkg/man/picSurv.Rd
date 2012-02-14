\name{picSurv}
\alias{picSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dataf, survObj, ncuts = 20, ...) 
{
    os <- order(dataf$size)
    os.surv <- (dataf$surv)[os]
    os.size <- (dataf$size)[os]
    psz <- tapply(os.size, as.numeric(cut(os.size, ncuts)), mean, 
        na.rm = TRUE)
    ps <- tapply(os.surv, as.numeric(cut(os.size, ncuts)), mean, 
        na.rm = TRUE)
    if (length(survObj@fit$coefficients) < 4) {
        plot(as.numeric(psz), as.numeric(ps), pch = 19, xlab = "Size at t", 
            ylab = "Survival to t+1", main = "Survival", ...)
        points(dataf$size[order(dataf$size)], surv(dataf$size[order(dataf$size)], 
            1, survObj), type = "l", col = 2)
    }
    else {
        plot(as.numeric(psz), as.numeric(ps), type = "n", pch = 19, 
            xlab = "Size at t", ylab = "Survival to t+1", main = "Survival", 
            ...)
        dataf$covariate <- as.factor(dataf$covariate)
        levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
        os.cov <- (dataf$covariate)[os]
        sizes <- dataf$size[!is.na(dataf$size)]
        sizes <- sizes[order(sizes)]
        ud <- unique(dataf$covariate)
        ud <- ud[!is.na(ud)]
        for (k in 1:length(ud)) {
            tp <- os.cov == ud[k]
            psz <- tapply(os.size[tp], as.numeric(cut(os.size[tp], 
                ncuts)), mean, na.rm = TRUE)
            ps <- tapply(os.surv[tp], as.numeric(cut(os.size[tp], 
                ncuts)), mean, na.rm = TRUE)
            points(as.numeric(psz), as.numeric(ps), pch = 19, 
                col = k)
            newd <- data.frame(size = sizes, size2 = sizes^2, 
                logsize = log(size), logsize2 = (log(size))^2, 
                covariate = rep(as.factor(ud[k]), length(sizes)))
            pred.surv <- predict(survObj@fit, newd, type = "response")
            points(newd$size, pred.surv, type = "l", col = k)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
