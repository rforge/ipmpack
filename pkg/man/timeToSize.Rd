\name{timeToSize}
\alias{timeToSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
timeToSize(starting.sizes, IPM, end.size, starting.env = 1, maxt = 100, prop.reach = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{starting.sizes}{
%%     ~~Describe \code{starting.sizes} here~~
}
  \item{IPM}{
%%     ~~Describe \code{IPM} here~~
}
  \item{end.size}{
%%     ~~Describe \code{end.size} here~~
}
  \item{starting.env}{
%%     ~~Describe \code{starting.env} here~~
}
  \item{maxt}{
%%     ~~Describe \code{maxt} here~~
}
  \item{prop.reach}{
%%     ~~Describe \code{prop.reach} here~~
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
function (starting.sizes, IPM, end.size, starting.env = 1, maxt = 100, 
    prop.reach = 0.01) 
{
    cutoff <- which(IPM@meshpoints > end.size, arr.ind = TRUE)[1]
    n.reach <- prop.reach * length(starting.sizes)
    if (IPM@n.env.class > 1) {
        if (length(starting.env) == 1) 
            starting.env <- rep(starting.env, length(starting.sizes))
        compound <- TRUE
        env.index <- IPM@env.index
        n.new.dist <- rep(0, length(IPM[1, ]))
        for (ev in 1:IPM@n.env.class) {
            index.new.dist <- findInterval(starting.sizes[starting.env == 
                ev], IPM@meshpoints) + 1
            index.new.dist[index.new.dist > length(IPM@meshpoints)] <- length(IPM@meshpoints)
            loc.sizes <- table(index.new.dist)
            n.new.dist[ev == IPM@env.index][as.numeric(names(loc.sizes))] <- loc.sizes
        }
        n.new.dist0 <- n.new.dist
    }
    else {
        compound <- FALSE
        index.new.dist <- findInterval(starting.sizes, IPM@meshpoints) + 
            1
        index.new.dist[index.new.dist > length(IPM@meshpoints)] <- length(IPM@meshpoints)
        loc.sizes <- table(index.new.dist)
        env.index <- rep(1, length(IPM@meshpoints))
        n.new.dist <- rep(0, length(IPM@meshpoints))
        n.new.dist[as.numeric(names(loc.sizes))] <- loc.sizes
        n.new.dist0 <- n.new.dist
    }
    ts.dist <- matrix(NA, length(n.new.dist), maxt)
    survivorship <- rep(NA, maxt)
    for (t in 1:maxt) {
        n.new.dist <- IPM@.Data \%*\% n.new.dist
        ts.dist[, t] <- n.new.dist
        if (!compound) {
            tot <- sum(n.new.dist[cutoff:length(IPM@meshpoints)])
            survivorship[t] <- sum(n.new.dist)/length(starting.sizes)
        }
        else {
            tot <- sumN <- 0
            for (ev in 1:IPM@n.env.class) {
                tot <- tot + sum(n.new.dist[env.index == ev][cutoff:length(IPM@meshpoints)])
                sumN <- sumN + sum(n.new.dist[env.index == ev])
            }
            survivorship[t] <- sumN/length(starting.sizes)
        }
        if (tot > n.reach) {
            time.reach <- t
            break()
        }
    }
    if (t == maxt) 
        time.reach <- maxt
    par(mfrow = c(2, 2), bty = "l")
    plot(IPM@meshpoints, n.new.dist0[env.index == 1], type = "l", 
        xlab = "size", ylab = "n", ylim = range(c(n.new.dist0, 
            n.new.dist)))
    points(IPM@meshpoints, n.new.dist[env.index == 1], type = "l", 
        col = 2)
    abline(v = IPM@meshpoints[cutoff], lty = 3)
    plot(survivorship[1:t], xlab = "Time", ylab = "Survivorship", 
        type = "l")
    if (time.reach > 5) {
        image(1:time.reach, IPM@meshpoints, t(log(ts.dist)), 
            xlab = "Time steps", ylab = "Number in each size class")
        contour(1:time.reach, IPM@meshpoints, t(log(ts.dist[, 
            1:time.reach])), add = TRUE)
    }
    print(paste("Time to reach:", time.reach))
    return(list(ts.dist = ts.dist, time.reach = time.reach, survivorship = survivorship))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
