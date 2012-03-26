\name{picSurv}
\alias{picSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Makes pictures of survival.
}
\description{
Produces figures of value for assessing survival fit given data.
}
\usage{
picSurv(dataf, survObj, ncuts = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and `sizeNext'(`size' is continuous stage variable at t, `sizeNext' is continuous stage at t+1); facultatively, dataf may include `covariate' and `covariateNext' for a single discrete covariate, indicating values at t, and at t+1, respectively; these must take values of sequential integers, starting at `1'. For models fitting growth increment, 'incr' or 'logincr' may be directly provided as a column in the dataframe, otherwise they are calculated as dataf$sizeNext-dataf$size or log(dataf$sizeNext-dataf$size), respectively.
  }
  \item{survObj}{
an object of class survObj.
}
  \item{ncuts}{
number of consecutive values for which means of survival and continuous (e.g. size) stage are taken for the plotting.
}
  \item{\dots}{
extra arguments to plot (e.g, ylim, etc).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns nothing.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeSurvObj}}, \code{\link{picGrow}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff <- generateData()
sv1 <- makeSurvObj(dff)
picSurv(dff,sv1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
