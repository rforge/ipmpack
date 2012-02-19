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
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean sizeNext directly; the density function around this is normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.declinevar")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean sizeNext directly; the density function around this is normal and the sd may change with size. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.Hossfeld")}}{
Methods to predict the probability density of \code{sizeNext} using the Hossfeld function. The prediction
supplies mean incr which is added to size to obtain sizeNext; the  density function  around this is normal. 
 }

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.incr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean incr which is added to size to obtain sizeNext; the density function around this is normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.incrdeclinevar")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean sizeNext directly; the density function around this is normal and the sd may change with size.
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov.incr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and covariates defined in data.frame. The prediction 
supplies mean incr which is added to size to obtain sizeNext; the density function around this is normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.logincr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean logincr which is exponentiated and added to size to obtain sizeNext; the density function around this is log normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj.logincrdeclinevar")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies logincr which is exponentiated and added to size to obtain sizeNext; the density function around this is lognormal; 
the sd may change with size. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov.logincr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and covariates defined in data.frame. The prediction 
supplies logincr which is exponentiated and added to size to obtain sizeNext; the density function around this is lognormal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and covariates defined in data.frame. The prediction 
supplies mean sizeNext directly; the density function around this is normal. 
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
