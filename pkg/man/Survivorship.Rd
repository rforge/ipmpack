\name{Survivorship}
\alias{Survivorship}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to estimate survivorship
}
\description{
  Calculates the fraction of the cohort surviving across age for a
  chosen starting size
}
\usage{
Survivorship(IPM.matrix, size1, maxage)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPM.matrix}{
    an IPM.matrix object describing growth and survival transitions
    across size and environment
  }  \item{size1}{
    a starting size for age 1
  }
  \item{maxage}{
    the maximum age up to which survivorships is desired for
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{surv.curv}{vector of length maxage providing survivorship at each
 age from 1 to maxage}
 \item{state.agesurv}{matrix of dimensions nBigMatrix*maxage providing
   the population structure at every age for a cohort starting with an
   individual of size size1
 }
 \item{mortality}{vector of length maxage providing mortality at each
 age from 1 to maxage}
}
\references{
Tuljapurkar & Horvitz, 2006. From stage to age in variable
environments. Life expectancy and survivorship. Ecology 87 p1497-1509
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LifeExpect}}
}
\examples{
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
su <- Survivorship(Tmatrix,1,300)
plot(su$surv.curv,type="l")


Tmatrix <-create.compound.Tmatrix(minSize=min(dff$size,na.rm=TRUE),
maxSize=max(dff$size,na.rm=TRUE),envMatrix=makeEnvObj(dff),
growObj=makeGrowthObj(dff, explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObj(dff,explanatoryVariables="size+size2+covariate"),DiscreteTrans=1)
su <- Survivorship(Tmatrix,1,300)
plot(su$surv.curv,type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
