\name{makeListTmatrix}
\alias{makeListTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Builds list of T matrices from list growth and survival objects.
}
\description{
  Uses posteriors (or other) lists of growth and survival objects to create a list of T matrices.
}
\usage{
makeListTmatrix(growObjList, survObjList, nBigMatrix, minSize, maxSize, cov = FALSE, envMat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growObjList}{
    list of growth objects.
  }
  \item{survObjList}{
    list of survival objects.
  }
  \item{nBigMatrix}{
    number of meshpoints.
  }
  \item{minSize}{
    minimum size. 
  }
  \item{maxSize}{
    maximum size.
  }
  \item{cov}{
    level of the covariate.
  }
  \item{envMat}{
    environmental matrix - defaults to NULL.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListTmatrix}}, \code{\link{makeListIPMs}}
}
\examples{
dff <- generateData()
svlist <- makePostSurvivalObjs(dataf = dff, meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
grlist <- makePostGrowthObjs(dataf = dff, nitt = 5000)

# Only first 10 for speed
Tlist <- makeListTmatrix(growObjList = grlist[1:10], survObjList = svlist[1:10], nBigMatrix = 20, minSize = -5, maxSize = 20, cov = FALSE, envMat = NULL)

# Plots four T matrices at random:
par(mfrow=c(2,2))
image(Tlist[[1]]@meshpoints, Tlist[[1]]@meshpoints, t(Tlist[[1]]), xlab = "Continuous (e.g. size) stage at time t", ylab = "Continuous (e.g. size) stage at time t+1")

image(Tlist[[2]]@meshpoints, Tlist[[2]]@meshpoints, t(Tlist[[2]]), xlab = "Continuous (e.g. size) stage at time t", ylab = "Continuous (e.g. size) stage at time t+1")

image(Tlist[[3]]@meshpoints, Tlist[[3]]@meshpoints, t(Tlist[[3]]), xlab = "Continuous (e.g. size) stage at t", ylab = "Continuous (e.g. size) stage at t+1")

image(Tlist[[4]]@meshpoints, Tlist[[4]]@meshpoints, t(Tlist[[4]]), xlab = "Continuous (e.g. size) stage at time t", ylab = "Continuous (e.g. size) stage at time t+1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
