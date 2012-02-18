\name{picGrow}
\alias{picGrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Function to make pictures of data with growth models
}
\description{
Takes the data-file and a growth object and shows the model fit over the data. 
}
\usage{
picGrow(dataf, growObj)
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

  \item{growObj}{
	an object of class growthObj that contains a fit for which R has methods for the function "predict"
}
}
\details{
Note that this model will ONLY work with growth objects that contain objects
of class glm or lm - i.e. ones for which R has defined methods for the function "predict"
}
\value{
Returns nothing
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeGrowthObjGeneral}}, \code{\link{picSurv}}
}
\examples{
dff<-generateData()
gr1 <- makeGrowthObjGeneral(dff)
picGrow(dff,gr1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
