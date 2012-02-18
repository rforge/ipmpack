\name{picSurv}
\alias{picSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
Function to make pictures of survival
}
\description{
Produces figures of value for assessing survival fit given data
}
\usage{
picSurv(dataf, survObj, ncuts = 20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `sizenext'(`size' is size at t, `sizenext' is size at t+1);
    facultatively, dataf may include `covariate'
    and `covariatel' for a single discrete covariate, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'. For models fitting growth increment, 'incr' or 'logincr' may be directly provided as
    a column in the dataframe, otherwise they are calculated as
    dataf$sizenext-dataf$size or log(dataf$sizenext-dataf$size),
  respectively. }
  \item{survObj}{
an object of class survObj
}
  \item{ncuts}{
number of consecutive values for which means of survival and size are taken for the plotting
}
  \item{\dots}{
extra arguments to plot (e.g, ylim, etc)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns nothing
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
\code{\link{makeSurvObjGeneral}}, \code{\link{picGrow}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff<-generateData()
sv1 <- makeSurvObjGeneral(dff)
picSurv(dff,sv1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
