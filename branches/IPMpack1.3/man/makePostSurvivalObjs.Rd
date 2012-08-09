\name{makePostSurvivalObjs}
\alias{makePostSurvivalObjs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make list of survival objects corresponding to posterior of Bayes fit
}
\description{
  Model takes a data-frame and fits a Bayesian model using the chosen explanatory variable, and the package MCMCglmm, fitting an error variance of 1; then create the resulting list of survival objects necessary to get Bayes posterior of IPM outputs
}
\usage{
makePostSurvivalObjs(dataf, explanatoryVariables="size+size2",
meanB = rep(0, 3), varB = rep(1e+10),burnin=3000, nitt = 50000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
    a dataframe with columns `size',`sizeNext', and `surv' (`size' is size at t, `sizeNext' is size at t+1, and `surv' is 0s and 1s, indicating if the individual
    survived or not); facultatively, dataf may include `covariate' (discrete or continuous) or other covariates. However, survival methods are not defined for a continuous covariate, or other covariates yet. 
  }
  \item{explanatoryVariables}{a formula describing the desired explanatory variables (interactions, etc) in classical R style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3' (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and `covariate'.}
  \item{meanB}{
    vector defining the means of the priors of the coefficients for survival (should be the same length as desired coeff). CURRENTLY NOT IMPLEMENTED.
  }
  \item{varB}{
    vector defining the variance of the prior of the coefficients for survival,as for the mean. CURRENTLY NOT IMPLEMENTED.
  }
  \item{burnin}{
    numeric defining burn in for MCMCglmm
}
  \item{nitt}{
    numeric defining the number of iterations
  }
}
\details{
Note that the use of priors is not yet implemented, and the function does not return output of Bayes analysis which might be useful for assessment; in development. Error variance is set to one. See Course Notes from MCMCpack for reasons that this results in slightly different parameters that would be obtained using a glm with family=binomial.
}
\value{
a list of survObjOverDisp survival objects, of length nitt; appropriate survival methods are defined that use Diggle et al. (2004)'s correction for estimating the expectation of survival following non-linear transform of a distribution on the logit scale.  
}
\references{
Diggle, P., P. Heagerty, K. Liang, and S. Zeger. 2004. Analysis of Longitudinal Data. 2nd edition. Oxford University Press.%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makePostGrowthObjs}}, \code{\link{makeListTmatrix}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff <- generateData()
svlist <- makePostSurvivalObjs(dataf=dff,
explanatoryVariables="size+size2+covariate",
meanB = rep(0, 3), varB = rep(1e+10), burnin = 100, nitt = 500)

plot(1:20,surv(1:20,data.frame(covariate=as.factor(1)),svlist[[1]]), type="n", xlab="size", ylab="survival")
for (k in 1:length(svlist)) 
  points(1:20,surv(1:20,data.frame(covariate=as.factor(1)),svlist[[k]]), type="l",col="grey")

}
