\name{surv-methods}
\docType{methods}
\alias{surv-methods}
\alias{surv,numeric,numeric,survObj-method}
\alias{surv,numeric,numeric,survObj.manycov-method}
\title{ ~~ Methods for Function surv  ~~}
\description{
 ~~ Methods for function \code{surv}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(size = "numeric", cov = "numeric", survObj = "survObj")}}{
Methods to predict probability of survival given a linear predictor based on 
various transforms of size and potentially a single discrete covariate. 
}

\item{\code{signature(size = "numeric", cov = "data.frame", survObj = "survObjMultiCov")}}{
Methods to predict probability of survival given a linear predictor based on 
various transforms of size and a dataframe of other covariates.  
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
