\name{makeListFmatrix}
\alias{makeListFmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build list of Fmatrices from lists of growth, survival and fertility objects
}
\description{
Uses posteriors (or other) lists of growth, fertility and survival objects to
create a list of Fmatrices
}
\usage{
makeListFmatrix(growObjList, survObjList, fecObjList, n.big.matrix, minsize, maxsize, cov = FALSE, env.mat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growObjList}{
list of growth objects
}
  \item{survObjList}{
list of survival objects
}
  \item{fecObjList}{
list of fertility objects
}
  \item{n.big.matrix}{
number of meshpoints
}
  \item{minsize}{
minimum size
}
  \item{maxsize}{
maximum size
}
  \item{cov}{
    level of covariate
  }
  \item{env.mat}{
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff<-generateData()
fvlist <- makePostFecObjs(dataf=dff, nitt=5000)
svlist <- makePostSurvivalObjs(dataf=dff, covsurv = FALSE,
meanB = rep(0, 3), varB = rep(1e+10), nitt = 5000)
grlist <- makePostGrowthObjs(dataf=dff, covgrow = FALSE,
meanB = rep(0, 3), varB = rep(1e+10), nitt = 5000)

#only first 10 for speed
Flist <- makeListFmatrix(growObjList=grlist[1:10],
survObjList=svlist[1:10],fecObjList=fvlist[1:10],
n.big.matrix=20, minsize=-5, maxsize=20, cov = FALSE, env.mat = NULL)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
