\name{create.MPM.Fmatrix}
\alias{create.MPM.Fmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to build discrete F matrix from the types of data structures used in IPMpack
}
\description{
A method that uses bins to construct a class discrete fertility transition matrix 
from the types of data structures used in IPMpack; assumes that there is a single fertility 
column in the data entitled fec that captures number of offspring produced by individuals over 
every size, and also assumes that all individuals are "born" into the first rows of the matrix 
defined by offspring.classes (where offspring.classes=1 means between bins[1] and bins[2], 
offspring.classes=2 means between bins[2] and bins[3], etc...) in proportions offspring.prop, 
which will be rescaled to sum to 1 if necessary  
}
\usage{
create.MPM.Fmatrix(dataf,bins, offspring.classes=c(1),offspring.prop,n.env=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size', `fec' (`size' is size at t, and `fec' indicates the number of
    seeds produced by that individual); facultatively, dataf may include `covariate'
    and `covariatel' for discrete covariates, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'. }
  \item{bins}{
    numeric, limits of the desired size classes
  }
  \item{n.env}{
    currently not used
  }
}
\details{
  The median number of offspring produced by individuals in each bin is
  identified. Fertility values are distributed across rows and columns accordingl. To
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
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.IPM.Fmatrix}}, \code{\link{create.MPM.Tmatrix}}
}
\examples{
dff<-generateData()
Fmatrix <-create.MPM.Fmatrix(dff,bins=c(1,4,10,15),offspring.classes=1, offspring.prop=1)
Fmatrix
Fmatrix <-create.MPM.Fmatrix(dff,bins=c(1,4,10,15),offspring.classes=c(1,2), offspring.prop=c(0.5,0.5))
Fmatrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
