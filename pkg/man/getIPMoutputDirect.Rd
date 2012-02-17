\name{getIPMoutputDirect}
\alias{getIPMoutputDirect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get IPM output from a list of growth, survival, and (facultatively)
fertility objects
}
\description{
Get summary statistics including life expectancy, passage time, and if a
fertility matrix is available, lambda, stable size, etc., for e.g.,
posterior distribution of survival and gowth objects, slightly swifter
approach than getIPMoutput since matrices are not stored 
}
\usage{
getIPMoutputDirect(survObjList, growObjList, target.size = c(),
n.big.matrix, minsize, maxsize, cov = FALSE, fecObjList = NULL,
env.mat = NULL, n.size.to.age = 0, size.start = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survObjList}{
list of survival objects
}
  \item{growObjList}{
list of growth objects
}
  \item{target.size}{
    Target size for passage time - if none is provided defaults to the
    median of the IPM meshpoints
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
    boolean identifying if a discrete covariate should be fitted
}
  \item{fecObjList}{
list of fertiltiy objects (facultative)
}
  \item{env.mat}{
    environmental matrix (facultative)
}
  \item{n.size.to.age}{
    numeric describing how many size to age defined (0 - 100s)
}
  \item{size.start}{
    starting size for size to age
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{LE}{matrix of life expectancies, columsn correspond to
  meshpoints, rows corresponding to each element of the list of Tmatrices}
  \item{ptime}{matrix of passage times to the target.size from each of
    the meshpoints (columns) and for each element in the Tmatrix list
    (columns)}
  \item{lambda}{vector of lambda corresonding to value obtained
    combining each element of the list of Tmatrices with the
    corresponding element in the list of Fmatrices; if no Fmatrix list
    is provided, this is a vector of NAs}
  \item{stable.size}{matrix of stable size distributions rows corresonding to values obtained
    combining each element of the list of Tmatrices with the
    corresponding element in the list of Fmatrices; if no Fmatrix list
    is provided, this is a matrix of NAs}
  \item{meshpoints} {matrix meshpoints}
  \item{resSize} {matrix providing target sizes for size to age estimate (assuming age=1
    at size.start), of length n.size.to.age space equally between the
    smallest and largest meshpoints}
  \item{resAge} {matrix providing time in time-steps to get to resSize,
    rows corresonding sequential elements in the list of Tmatrices} 
  \item{surv.par} {matrix of posteriors of survival parameters (columns,
    each sample is a row)}
  \item{grow.par} {matrix of posteriors of growth parameters (in columns)} 
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
dff <- generateData()

svlist <- makePostSurvivalObjs(dataf=dff, 
meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
grlist <- makePostGrowthObjs(dataf=dff,  nitt=5000)
fvlist <- makePostFecObjs(dataf=dff, nitt=5000)

rc <- getIPMoutputDirect(survObjList=svlist[1:10], growObjList=grlist[1:10], target.size = 10,
n.big.matrix=20, minsize=-2, maxsize=15, cov = FALSE, fecObjList =fvlist[1:10],
env.mat = NULL, n.size.to.age = 5, size.start = -1)

names(rc)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
