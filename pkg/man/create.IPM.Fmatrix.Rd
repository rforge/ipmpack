\name{create.IPM.Fmatrix}
\alias{create.IPM.Fmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to build F matrices
}
\description{
Uses fertility objects to construct a matrix defining probabailities for
transitions between stages due to fecundity
}
\usage{
create.IPM.Fmatrix(n.env.class = 1, n.big.matrix = 50,
minsize = -1, maxsize = 50, chosen.cov = 1,  fecObj, integrate.type,correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.env.class}{
    numeric, number of environmental classes, always = 1 for
    non-compound matrices
  }
  \item{n.big.matrix}{
    numeric, number of size bins in the T matrix, defaults to 50
  }
  \item{minsize}{
    numeric, minimum size of the T matrix, defaults to -1
  }
  \item{maxsize}{
    numeric, maximum size of the T matrix, defaults to -1
  }
  \item{chosen.cov}{
    numeric indicating level of the discrete covariate, or vector
    setting values where multiple covarites are modeled
  }
  \item{fecObj}{
    fertility object
  }
  \item{integrate.type}{
    integration type, defaults to "midpoint" (which uses prob. density
    function); other option is "cumul" (which uses the cumul. density function)
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will
    multiply every column of the IPM by a constant sufficient to adjust
    values to those predicted for survival at that size. 
  }
 
}
\details{
%%
}
\value{
  an object of class IPM.matrix of dimensions n.big.matrix or n.big.matrix+ndiscrete classes (defined by fecObj@offspring.splitter-1)
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
\code{\link{create.IPM.Tmatrix}}
}
\examples{


dff<-generateData()
fv1 <- makeFecObjGeneral(dff,
explanatoryVariables="size",Transform="log")
Fmatrix <-create.IPM.Fmatrix(n.big.matrix=20,minsize=min(dff$size, na.rm=TRUE),
maxsize=max(dff$size, na.rm=TRUE),fecObj=fv1)

slotNames(Fmatrix)

image(Fmatrix@meshpoints,Fmatrix@meshpoints,Fmatrix, xlab="Size at t+1",
ylab="Size at t")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
