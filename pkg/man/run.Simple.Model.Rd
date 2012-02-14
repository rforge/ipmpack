\name{run.Simple.Model}
\alias{run.Simple.Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to run a simple model
}
\description{
This function takes a dataframe, a chosen size, a min and max, and a
number of bins, and creates a Tmatrix and plots passage time and life expectancy.
}
\usage{
run.Simple.Model(dataf, chosen.size, minsize = c(), maxsize = c(), n.big.matrix = 500, do.plot = TRUE, is.log = TRUE,  integrate.type = "midpoint")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
  a data-frame with columns size, sizel (which is size in the NEXT time
  step), and surv, at a minimum%%     ~~Describe \code{dataf} here~~
}
  \item{chosen.size}{
the target size
    %%     ~~Describe \code{chosen.size} here~~
}
  \item{minsize}{
%%     ~~Describe \code{minsize} here~~
}
  \item{maxsize}{
%%     ~~Describe \code{maxsize} here~~
}
  \item{n.big.matrix}{
%%     ~~Describe \code{n.big.matrix} here~~
}
  \item{do.plot}{
%%     ~~Describe \code{do.plot} here~~
}
  \item{is.log}{
%%     ~~Describe \code{is.log} here~~
}

  \item{integrate.type}{
%%     ~~Describe \code{integrate.type} here~~
}
}
\details{
This function is a utility for running a whole bunch of species... or
something... %%  ~~ If necessary, more details than the description above ~~
}
\value{
  a list including LE which is life
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
dff<-generateData()
run.Simple.Model(dataf=dff, chosen.size=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
