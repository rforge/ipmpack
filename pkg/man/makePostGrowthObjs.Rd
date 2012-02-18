\name{makePostGrowthObjs}
\alias{makePostGrowthObjs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make list of growth objects corresponding to posterior of Bayes fit
}
\description{
  Model takes a data-frame and fits a Bayesian model to the chosen
  "responseType" (sizenext, incr, etc) using the chosen explanatory
  variable, and the package MCMCglmm; then
  create the resulting list of growth objects necessary to get Bayes
  posterior of IPM outputs
}
\usage{
makePostGrowthObjs(dataf, explanatoryVariables,responseType,
meanB = rep(0, 3), varB = rep(1e+10), nitt = 50000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size',
    and `sizel' (`size' is size at t, `sizel' is size at t+1); facultatively, dataf may include `covariate'
    (discrete or continuous) or other covariates
  }
  \item{explanatoryVariables}{a formula describing the desired
    explanatory variables (interactions, etc) in classical R
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate'.}
  \item{responseType}{possible values include `sizenext' (size in one
    discrete time step) and `incr' and `logincr'}  
  \item{meanB}{
    vector defining the means of the priors of the coefficients for growth (should be the same length as desired coeff)
  }
  \item{varB}{
    vector defining the variance of the prior of the coeff for growth,
    as for the mean
  }
  \item{nitt}{
    numeric defining the number of iterations
}

}
\details{
}
\value{
  a list of growth objects
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
\code{\link{makePostSurvivalObjs}}, makePostSurvObjs 
}
\examples{
dff<-generateData()
grlist <- makePostGrowthObjs(dataf=dff,
			explanatoryVariables="size+size2+covariate",
            responseType="sizenext",
            meanB=rep(0,3),varB=rep(1e10),nitt = 5000)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
