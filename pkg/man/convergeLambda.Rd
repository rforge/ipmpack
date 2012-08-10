\name{convergeLambda}
\alias{convergeLambda}
\alias{convergeR0}
\alias{convergeLifeExpectancyFirstBin}
\alias{convergeLifeExpectancyLastBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Iterates until obtaining the number of bins required so that the difference in the chosen measure (lambda, R0, life expectancy of the first bin) falls below a chosen tolerance level
}
\description{
  Increases bin number by a specified \code{binIncrease} until the difference in the chosen measure falls below a tolerance level 
}
\usage{
convergeLambda(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
			discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-4,binIncrease=5)
convergeR0(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
			discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-4,binIncrease=5)
convergeLifeExpectancyFirstBin(growObj, survObj,nBigMatrix, minSize, maxSize, 
				discreteTrans = 1, integrateType = "midpoint", correction = "none", 
				tol=1e-1,binIncrease=5)
convergeLifeExpectancyLastBin(growObj, survObj,nBigMatrix, minSize, maxSize, 
				discreteTrans = 1, integrateType = "midpoint", correction = "none", 
				tol=1e-1,binIncrease=5)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growObj}{
    a growth object.
  }
  \item{survObj}{
    a survival object.
  }
  \item{fecObj}{
    a fecundity object.
  }
  \item{nBigMatrix}{
    numeric, initial number of bins of size used in the matrix - will be increased for the assessment
  }
  \item{minSize}{
    numeric, minimum size used for meshpoints
  }
  \item{maxSize}{
    numeric, maximum size used for meshpoints of the T matrix.
  }
  \item{discreteTrans}{
	matrix of discrete transitions; or 1 if there is none
  }
  \item{integrateType}{
    integration type.
  }
  \item{correction}{
    correction (see createIPMTmatrix)
  }
  \item{preCensus}{
    boolean defining whether fecundity is pre or post census; defaults to pre
  }
  \item{tol}{
   desired tolerance level 
  }
  \item{binIncrease}{
   increments in increase in the number of bins (should be an integer
  }

  
}
\value{
  \item{binIncrease}{the number of bins used to increase matrix size in assessing tolerance}
  \item{Tmatrix}{the final Tmatrix if only LE is being considered}
  \item{IPM}{the final IPM}
  \item{R0}{the final R0}
  \item{lambda}{the final lambda}
  \item{LE}{the final vector of life expectancies}
}
\details{
For convergeLifeExpectancyFirstBin, if discrete stages are included via discreteTrans then convergeLifeExpectancyFirstBin will use the first discrete bin. 
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
This code was modified from original code by Melissa Eitzel. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diagnosticsTmatrix}}
}
\examples{
dff<-generateData()
gr1<-makeGrowthObj(dff)
sv1<-makeSurvObj(dff)
fv1<-makeFecObj(dff,Transform="log")

res <- convergeLambda(growObj=gr1, survObj=sv1, fecObj=fv1, nBigMatrix=10, minSize=-2, maxSize=15, 
			discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-3,binIncrease=10)


res <- convergeR0(growObj=gr1, survObj=sv1, fecObj=fv1, nBigMatrix=10, minSize=-2, maxSize=15, 
			discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-3,binIncrease=10)

res <- convergeLifeExpectancyFirstBin(growObj=gr1, survObj=sv1, nBigMatrix=10, minSize=-2, maxSize=15, 
				discreteTrans = 1, integrateType = "midpoint", correction = "none",  tol=1e-1,binIncrease=10)

res <- convergeLifeExpectancyLastBin(growObj=gr1, survObj=sv1, nBigMatrix=10, minSize=-2, maxSize=15, 
				discreteTrans = 1, integrateType = "midpoint", correction = "none",  tol=1e-1,binIncrease=10)

}
