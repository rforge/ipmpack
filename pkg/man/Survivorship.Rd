\name{Survivorship}
\alias{Survivorship}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimates survivorship between two time censuses.
}
\description{
  Calculates the fraction of the cohort surviving across age for a chosen starting continuous stage value.
}
\usage{
Survivorship(IPM.matrix, size1, maxAge)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPM.matrix}{
    an IPM.matrix object describing growth and survival transitions across stage (e.g. size) and environment.
  }  \item{size1}{
    a starting size for age 1.
  }
  \item{maxAge}{
    the maximum age up to which survivorship is desired for or possible.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{surv.curv}{vector of length maxAge providing survivorship at each
 age from 1 to maxAge.}
 \item{stageAgeSurv}{matrix of dimensions nBigMatrix*maxAge providing
   the population structure at every age for a cohort starting with an
   individual of size size1.
 }
 \item{mortality}{vector of length maxAge providing mortality at each
 age from 1 to maxAge.}
}
\references{
Tuljapurkar & Horvitz. 2006. From stage to age in variable
environments. Life expectancy and survivorship. Ecology 87, p1497-1509.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LifeExpect}}
}
\examples{
# For only continuous stages (e.g. size)
dff <- generateData()
Tmatrix <- create.IPM.Tmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff), survObj = makeSurvObj(dff))
su <- Survivorship(Tmatrix, 1, 300)
plot(su$surv.curv, type = "l", col = "dark gray", ylab = "Survivorship", xlab= "Continuous (e.g. size) stage", ylim = c(0,1))

# For continuous (e.g. size) and discrete (e.g. seedbank) stages
Tmatrix <- create.compound.Tmatrix(minSize = min(dff$size,na.rm = TRUE),
maxSize = max(dff$size,na.rm = TRUE), envMatrix = makeEnvObj(dff),
growObj = makeGrowthObj(dff, explanatoryVariables = "size+size2+covariate"),
survObj = makeSurvObj(dff, explanatoryVariables = "size+size2+covariate"), discreteTrans = 1)
su <- Survivorship(Tmatrix,1,300)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
