\name{makePostSurvivalObjs}
\alias{makePostSurvivalObjs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make list of survival objects corresponding to posterior of Bayes fit
}
\description{
  Model takes a data-frame and fits a Bayesian model using the chosen explanatory
  variable, and the package MCMCglmm; then create the resulting list of
  survival objects necessary to get Bayes posterior of IPM outputs
}
\usage{
makePostSurvivalObjs(dataf, explanatoryVariables,
meanB = rep(0, 3), varB = rep(1e+10), nitt = 50000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size',
    `sizenext', and `surv' (`size' is size at t, `sizenext' is size at t+1,
    and `surv' is 0s and 1s, indicating if the individual
    survived or not); facultatively, dataf may include `covariate'
    (discrete or continuous) or other covariates
  }
  \item{explanatoryVariables}{a formula describing the desired
    explanatory variables (interactions, etc) in classical R
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate'.}
  \item{meanB}{
    vector defining the means of the priors of the coefficients for survival (should be the same length as desired coeff)
  }
  \item{varB}{
    vector defining the variance of the prior of the coeff for survival,
    as for the mean
  }
  \item{nitt}{
    numeric defining the number of iterations
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a list of survival objects, of length nitt
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
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
svlist <- makePostSurvivalObjs(dataf=dff,
explanatoryVariables="size+size2+covariate",
meanB = rep(0, 3), varB = rep(1e+10), nitt = 5000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
