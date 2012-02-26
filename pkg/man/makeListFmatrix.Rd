\name{makeListFmatrix}
\alias{makeListFmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build list of F (sexual reproduction) matrices from lists of growth, survival and fecundity objects
}
\description{
Uses posteriors (or other) lists of growth, fecundity and survival objects to
create a list of Fmatrices
}
\usage{
makeListFmatrix(fecObjList, nBigMatrix, minSize, maxSize, cov = FALSE, envMat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{fecObjList}{
	list of fecundity objects
}
  \item{nBigMatrix}{
	number of meshpoints
}
  \item{minSize}{
	minimum size
}
  \item{maxSize}{
	maximum size
}
  \item{cov}{
    level of covariate
  }
  \item{envMat}{
	environmental matrix, defaults to NULL
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list of Fmatrices. 
}
\references{
Yang, Jongejans, Yang & Bishop. 2011. The effect of consumers and mutualists of Vaccinum membranaceum at Mount St. Helens: dependence on successional context. PLoS One 10 p1-11

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeListTmatrix}}, \code{\link{makeListIPMs}}
}
\examples{
dff<-generateData()
fvlist <- makePostFecObjs(dataf=dff, nitt=5000)

#only first 10 for speed
Flist <- makeListFmatrix(fecObjList=fvlist[1:10],
nBigMatrix=20, minSize=-5, maxSize=20, cov = FALSE, envMat = NULL)

#make pics of two of them at random...
image(Flist[[1]]@meshpoints, Flist[[1]]@meshpoints, t(Flist[[1]]), 
xlab="size", ylab="size")

image(Flist[[3]]@meshpoints, Flist[[3]]@meshpoints, t(Flist[[3]]), 
xlab="size", ylab="size")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
