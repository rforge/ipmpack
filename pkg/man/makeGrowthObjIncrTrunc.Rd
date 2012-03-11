\name{makeGrowthObjIncrTrunc}
\alias{makeGrowthObjIncrTrunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to build growth objects for constant environments, or
  discretely varying environments, where incr follows a truncated normal distribution.
}
\description{
 A function that fits regressions defining growth (following next size,
  size increment, or log size increment) and from these build growth objects for which methods to build an IPM
  object are defined.
  }
\usage{
makeGrowthObjIncrTrunc(dataf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `sizeNext'(`size' is size at t, `sizeNext' is size at t+1);
    facultatively, dataf may include `covariate'
    and `covariatel' for a single discrete covariate, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'. For models fitting growth increment, 'incr' or 'logincr' may be directly provided as
    a column in the dataframe, otherwise they are calculated as
    dataf$sizeNext-dataf$size or log(dataf$sizeNext-dataf$size),
  respectively. }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class growthObjTruncIncr, an S4 object which contains the slots: 
   \item{fit}{an object of class \code{lm} or \code{glm} or \code{gls} that can be used with predict in the survival methods}
Slots can be listed by using slotNames(growthObj)   
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
Uses the package cenReg
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#generate data
dff<-generateData()
#make simple linear regression growth object relating size to size at t+1
#(gr1 <- makeGrowthObjIncrTrunc(dataf=dff))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
