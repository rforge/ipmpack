\name{create.IPM.Tmatrix}
\alias{create.IPM.Tmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to build T matrices
}
\description{
  Uses growth and survival objects to construct a matrix
  defining probabailities for transitions between sizes due to growth
  and survival. Extensions for transition to discrete classes are possible
}
\usage{
create.IPM.Tmatrix(n.env.class = 1, n.big.matrix = 50,
minsize = -1, maxsize = 50, chosen.cov = 1, growObj, survObj, integrate.type, correction)
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
    numeric, maximum size of the T matrix, defaults to 50
  }
  \item{chosen.cov}{
    numeric indicating level of the discrete covariate, or vector
    setting values where multiple covarites are modeled
  }
  \item{growObj}{
    growth object
  }
  \item{survObj}{
    survival object
  }
  \item{DiscreteTrans}{
    object of class DiscreteTrans, or numeric
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
  The number of bins (n.big.matrix) is combined with the minimum and
  maximum size to define the meshpoints of the IPM. Bins should be
  sufficient and the size range should encompass the size range of the
  data. If a "DiscreteTrans" exists, then discrete stages will be added
  to the IPM structure. If multiple discrete covariate levels are available,
  chosen.cov identifies the covariate value for which an IPM is
  required; if a series of covariates are being modeled, chosen.cov is a
  vector of these covariates, and growth will reflect these values.  
}
\value{
  an object of class IPM.matrix with dimensions n.big.matrix*n.big.matrix, or
  nbig.matrix+nrow(DiscreteTrans)
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
\code{\link{create.IPM.Fmatrix}}, \code{\link{diagnosticsTmatrix}}, \code{\link{makeDiscreteTrans}}  
}
\examples{
dff<-generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObjGeneral(dff), survObj=makeSurvObjGeneral(dff))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
