\name{create.MPM.Fmatrix}
\alias{create.MPM.Fmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to build discrete F matrix
}
\description{
A method that uses bins to construct a class discrete
fertility transition matrix 
}
\usage{
create.MPM.Fmatrix(dataf,bins,p.est=1, n.env=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size',
    `age', `fec' (`size' is size at t, `age' indicates if the
    individual is in its first year, and `fec' indicates the number of
    seeds produced by that individual); facultatively, dataf may include `covariate'
    and `covariatel' for discrete covariates, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'. }
  \item{bins}{
    numeric, limits of the desired size classes
  }
  \item{p.est}{
    numeric, probability of seed establishment
  }
  \item{n.env}{
    currently not used
  }
}
\details{
  The median number of seeds produced by individuals in each bin is
  identified, and the proportion of one year olds in each bin is
  determined. Fertility values are distributed across rows and columns accordingly; p.est allows
  adjustment of the fraction of seeds establishing. To
  get a classic Matrix Population model this can be added to values
  defined by create.MPM.Tmatrix.    
}
\value{
  an object of class IPM.matrix
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
create.IPM.Fmatrix, create.MPM.Fmatrix
}
\examples{
dff<-generateData()
Fmatrix <-create.MPM.Fmatrix(dff,bins=c(1,4,10,15))
Fmatrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
