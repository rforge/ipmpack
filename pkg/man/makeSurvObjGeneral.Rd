\name{makeSurvObjGeneral.Rd}
\alias{makeSurvObjGeneral.Rd}
\alias{picSurv}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to build survival objects}
\description{
  A function to fit logistic regressions defining
  survival following user defined forulas (e.g., size+size^2, etc)
  to build survival objects for which methods to build an IPM  object are defined.
}
\usage{
makeSurvObjGeneral(dataf,explanatoryVariables="size+size2+size3")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `surv'(`size' is size at t, `surv' is 0 for death of the individual
    and 1 for survival);
    facultatively, dataf may include `covariate'
    and `covariatel' for a single discrete covariate, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'.  }
  \item{explanatoryVariables}{a formula describing the desired
    explanatory variables (interactions, etc) in classical R
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate'.}
}

\details{

}
\value{
  An object of class survObj
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
  See manual for details on building case-specific survival objects. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{makeSurvObjGeneralManyCov}}
}
\examples{
#generate data
dff<-generateData()
#make simple logistic regression survival object relating survival to size at t
sv1 <- makeSurvObjGeneral(dff, explanatoryVariables="size")
#assess fit for model with discrete environmental classes fitted
sv1 <- makeSurvObjGeneral(dff, explanatoryVariables="size+covariate")
picSurv(dff,sv1,ncuts=30)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
