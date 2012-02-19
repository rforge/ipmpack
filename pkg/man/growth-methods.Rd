\name{growth-methods}
\docType{methods}
\alias{growth-methods}
\alias{growth,numeric,numeric,numeric,growthObj-method}
\alias{growth,numeric,numeric,numeric,growthObj.declinevar-method}
\alias{growth,numeric,numeric,numeric,growthObj.Hossfeld-method}
\alias{growth,numeric,numeric,numeric,growthObj.incr-method}
\alias{growth,numeric,numeric,numeric,growthObj.incrdeclinevar-method}
\alias{growth,numeric,numeric,numeric,growthObj.incrmanycov-method}
\alias{growth,numeric,numeric,numeric,growthObj.logincr-method}
\alias{growth,numeric,numeric,numeric,growthObj.logincrdeclinevar-method}
\alias{growth,numeric,numeric,numeric,growthObj.logincrmanycov-method}
\alias{growth,numeric,numeric,numeric,growthObj.manycov-method}
\title{ Methods for Function growth  ~~}
\description{
 Methods for defining the growth component of the IPM
}
\section{Methods}{
\describe{

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj")}}{
Methods to predict the probability density of \code{sizeNext} using a range of transforms of the current size, 
and potentially one discrete covariate in the predictor, which supplies mean sizeNext directly; the 
density function is normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.declinevar")}}{
Methods to predict the probability density of \code{sizeNext} using a range of transforms of the current size, 
and potentially one discrete covariate in the predictor, which supplies mean sizeNext directly; the 
density function is normal, but the sd may change with size. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.Hossfeld")}}{
Methods to predict the probability density of \code{sizeNext} using the Hossfeld function
 which supplies mean incr which is added to size to obtain sizeNext; the  density function is normal. 
 }

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.incr")}}{
Methods to predict the probability density of \code{sizeNext} using a range of transforms of the current size, 
and potentially one discrete covariate in the predictor, which supplies mean incr which is added to size to 
obtain sizeNext; the density function is normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.incrdeclinevar")}}{
Methods to predict the probability density of \code{sizeNext} using a range of transforms of the current size, 
and potentially one discrete covariate in the predictor, which supplies mean incr which is added to size to 
obtain sizeNext; the density function is normal, but the sd may change with size. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov.incr")}}{
Methods to predict the probability density of \code{sizeNext} using a range of transforms of the current size, 
and a data-frame of covariates in the predictor, which supplies mean incr which is added to size to 
obtain sizeNext; the density function is normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.logincr")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.logincrdeclinevar")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjMultiCov.logincr")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjMultiCov")}}{
%%  ~~describe this method here~~
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
