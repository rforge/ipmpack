\name{makeListTmatrix}
\alias{makeListTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Build list of Tmatrices from list growth and survival objects
}
\description{
  Uses posteriors (or other) lists of growth and survival objects to
  create a list of Tmatrices
}
\usage{
makeListTmatrix(growObjList, survObjList, n.big.matrix, minsize, maxsize, cov = FALSE, env.mat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growObjList}{
    list of growth objects
  }
  \item{survObjList}{
    list of survival objects
  }
  \item{n.big.matrix}{
    number of meshpoints
  }
  \item{minsize}{
    min size 
  }
  \item{maxsize}{
    max size
  }
  \item{cov}{
    level of the covariate
  }
  \item{env.mat}{
    environmental matrix - defaults to NULL
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
%dff<-generateData()
%svlist <- makePostSurvivalObjs(dataf=dff, covsurv = FALSE,
%meanB = rep(0, 3), varB = rep(1e+10), nitt = 1000)
%grlist <- makePostGrowthObjs(dataf=dff, covgrow = FALSE, nitt=5000)

#use only first 10 for speed
%Tlist <- makeListTmatrix(growObjList=grlist[1:10], survObjList=svlist[1:10],
%n.big.matrix=20, minsize=-5, maxsize=20, cov = FALSE, env.mat = NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
