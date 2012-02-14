\name{convertIncrement}
\alias{convertIncrement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert size increment according to time elapsed
}
\description{
In many data-sets the time-interval between sequential measurements will
be slightly different; this provides a simple method to adjust them. 
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
%%  ~~  none If necessary, more details than the description above ~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
