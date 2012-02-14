\name{plotResultsStochStruct}
\alias{plotResultsStochStruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotResultsStochStruct(tvals, st, covtest, n.runin = 15, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tvals}{
%%     ~~Describe \code{tvals} here~~
}
  \item{st}{
%%     ~~Describe \code{st} here~~
}
  \item{covtest}{
%%     ~~Describe \code{covtest} here~~
}
  \item{n.runin}{
%%     ~~Describe \code{n.runin} here~~
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
function (tvals, st, covtest, n.runin = 15, ...) 
{
    if (length(st$IPM.here@meshpoints) < length(st$rc[, 1])) {
        do.seed <- TRUE
        index.notseed <- (1:length(st$IPM.here@meshpoints)) + 
            1
    }
    else {
        do.seed <- FALSE
        index.notseed <- (1:length(st$IPM.here@meshpoints))
    }
    par(mfrow = c(2, 2), bty = "l")
    plot(tvals[n.runin:length(tvals)], colSums(st$rc[index.notseed, 
        n.runin:length(tvals)]), log = "y", xlab = "Time (years)", 
        ylab = "Population size (plants)", type = "l", ...)
    abline(v = 1:max(tvals), lty = 3)
    covtestplot <- exp(mean(log(colSums(st$rc[index.notseed, 
        n.runin:length(tvals)]))) + ((covtest - mean(covtest))/sd(covtest)) * 
        sd(log(colSums(st$rc[index.notseed, n.runin:length(tvals)]))))
    points(tvals, covtestplot + 1, type = "l", lty = 3, col = 2)
    if (do.seed) {
        plot(tvals[n.runin:length(tvals)], st$rc[1, n.runin:length(tvals)], 
            log = "y", xlab = "Time (years)", ylab = "Population size (seeds)", 
            type = "l", ...)
        abline(v = 1:max(tvals), lty = 3)
        covtestplot <- exp(mean(log(st$rc[1, n.runin:length(tvals)])) + 
            ((covtest - mean(covtest))/sd(covtest)) * sd(log(st$rc[1, 
                n.runin:length(tvals)])))
        points(tvals, covtestplot + 1, type = "l", lty = 3, col = 2)
    }
    image(tvals[n.runin:length(tvals)], st$IPM.here@meshpoints, 
        t(log(st$rc[index.notseed, n.runin:length(tvals)] + 1)), 
        ylab = "Size (e.g. leaf length)", xlab = "Time", ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
