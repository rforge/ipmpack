\name{makePostGrowthObjs}
\alias{makePostGrowthObjs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make list of growth objects corresponding to posterior of Bayes fit
}
\description{
  Model takes a data-frame and fits a Bayesian model to the chosen
  "responseType" (sizeNext, incr, etc) using the chosen explanatory
  variable, and the package MCMCglmm; then creates the resulting list of growth objects necessary to get Bayes
  posteriors of IPM outputs. Note that priors are not yet implemented and growth methods for multiple covariates likewise.  
 }
\usage{
makePostGrowthObjs(dataf, explanatoryVariables,responseType,
meanB = rep(0, 3), varB = rep(1e+10), burnin=3000,nitt = 50000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size',
    and `sizeNext' (`size' is size at t, `sizeNext' is size at t+1); facultatively, dataf may include `covariate'
    (discrete or continuous) or other covariates
  }
  \item{explanatoryVariables}{a formula describing the desired
    explanatory variables (interactions, etc) in classical R
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate'.}
  \item{responseType}{possible values include `sizeNext' (size in one
    discrete time step) and `incr' and `logincr'}  
  \item{meanB}{
    vector defining the means of the priors of the coefficients for growth (should be the same length as desired coeff). CURRENTLY NOT IMPLEMENTED.
  }
  \item{varB}{
    vector defining the variance of the prior of the coeff for growth,
    as for the mean. CURRENTLY NOT IMPLEMENTED.
  }
  \item{burnin}{
    numeric defining burn in for MCMCglmm
}
  \item{nitt}{
    numeric defining the number of iterations
}

}
\details{
Note that the use of priors is not yet implemented, and the function does not return output
of Bayes analysis which might be useful for assessment; in development. 

}
\value{
  a list of growth objects
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makePostSurvivalObjs}},  \code{\link{makeListPmatrix}}
}
\examples{
dff <- generateData()
grlist <- makePostGrowthObjs(dataf=dff,
			explanatoryVariables="size+size2+covariate",
            responseType="sizeNext",
            meanB=rep(0,3),varB=rep(1e10),burnin=100,nitt = 500)
            
grlist[[1]]
}

