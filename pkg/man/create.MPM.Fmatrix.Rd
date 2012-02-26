\name{create.MPM.Fmatrix}
\alias{create.MPM.Fmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds discrete F matrices from the types of data structures used in IPMpack.
}
\description{
A method that uses bins to construct a class discrete fecundity matrix population model (F) from the types of data structures used in IPMpack. It assumes that there is a single fecundity column in the data entitled "fec" that captures number of offspring produced by individuals over every continuous stage (e.g. size). Individuals can be "born" into different stages of the matrix (e.g. seedbank, seedling) defined by offspringClasses (where offspringClasses=1 means between bins[1] and bins[2], offspringClasses=2 means between bins[2] and bins[3], etc...) in proportions offspringProp, which will be rescaled to sum to 1 if necessary. 
}
\usage{
create.MPM.Fmatrix(dataf,bins, offspringClasses=1, offspringProp=1, nEnv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    A dataframe with columns `size', `fec' (`size' is size at t, and `fec' indicates the number of propagules produced by that individual); facultatively, dataf may include `covariate' and `covariateNext' for discrete covariates, indicating values at time t, and at t+1, respectively. These covariates must take values of sequential integers, starting at `1'.
    }
  \item{bins}{
    numeric, limits of the desired discrete (size) classes.
  }
  \item{nEnv}{
    currently not used.
  }
  \item{offspringClasses}{
  the bins of the matrix in which offspring appear.
  }
  \item{offspringProp}{
  the proportion of offspring appearing in each bin.
  }
}
\details{
  The median number of offspring produced by individuals in each bin is
  identified. Fecundity values are distributed across rows and columns accordingly. To
  get a classic matrix population model (A) this can be added to values defined by create.MPM.Tmatrix (A=F+T).    
}
\value{
  an object of class IPM.matrix.
}
\references{
For information on F matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.

For information on collapsing large matrices: Salguero-Gomez & Plotkin. 2010. Matrix dimensions bias demographic inferences: implications for comparative plant demography. The American Naturalist 176, p710-720.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.IPM.Fmatrix}}, \code{\link{create.MPM.Tmatrix}}
}
\examples{
dff <- generateData()
matrixCuts <- c(1, 4, 10, 15)
Fmatrix <-create.MPM.Fmatrix(dff, bins=matrixCuts, offspringClasses = c(1,2), offspringProp = c(0.5,0.5))
Fmatrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
