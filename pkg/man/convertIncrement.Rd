\name{convertIncrement}
\alias{convertIncrement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert size increment according to time elapsed
}
\description{
This function adjust the intervals of census data that is not annual. 
}
\usage{
convertIncrement(dataf, nyrs = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a data-frame with columns `size' `sizel`, `exactDate', `exactDatel' 
}
  \item{nyrs}{
    the number of years between sequential measurements (i.e. if census
    interval is 5, nyrs=5)
}
}
\details{
In many data-sets the time-interval between census measurements will be greater than  a single year and/or vary across intervals. This function takes a data frame \code{dataf} and uses columns with the term \code{extractData} to adjust the size increment to the number of years given in \code{nyears}.  Defaults to annual.  %%  ~~  none  If necessary, more details than the description above ~~
}
\value{
Returns the adjusted increments.%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
