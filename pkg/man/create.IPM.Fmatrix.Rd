\name{create.IPM.Fmatrix}
\alias{create.IPM.Fmatrix}
%- Also NEED an '\alias' for EACH other topic documented here. 
\title{
  Function to build F matrices
}
\description{
Uses fecundity objects to construct a matrix defining per-capita contribution to recruitment stages (e.g., propagules [seed, spore], seedlings, calfs) by reproductive stages due to sexual reproduction
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
    setting values where multiple covariates are modeled
  }
  \item{fecObj}{
    fecundity object
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
For information on F matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112

For midpoint: Zuidema, Jongejans, Chien, Duirng & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98 p345-355

For multiple-vital rate integration on fecundity: Yang, Jongejans, Yang & Bishop. 2011. The effect of consumers and mutualists of Vaccinum membranaceum at Mount St. Helens: dependence on successional context. PLoS One 10 p1-11

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
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
fv1 <- makeFecObj(dff,explanatoryVariables="size",Transform="log")
Fmatrix <-create.IPM.Fmatrix(n.big.matrix=20,minsize=min(dff$size, na.rm=TRUE),maxsize=max(dff$size, na.rm=TRUE),fecObj=fv1)

slotNames(Fmatrix)

image(Fmatrix@meshpoints,Fmatrix@meshpoints,t(Fmatrix), xlab="Size at t",ylab="Size at t+1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
