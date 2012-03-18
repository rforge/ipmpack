\name{growth-methods}
\docType{methods}
\alias{growth-methods}
\alias{growth,numeric,numeric,data.frame,growthObjMultiCov-method}
\alias{growth,numeric,numeric,data.frame,growthObjMultiCov.incr-method}
\alias{growth,numeric,numeric,data.frame,growthObjMultiCov.logincr-method}
\alias{growth,numeric,numeric,numeric,growthObj-method}
\alias{growth,numeric,numeric,numeric,growthObjDeclineVar-method}
\alias{growth,numeric,numeric,numeric,growthObjHossfeld-method}
\alias{growth,numeric,numeric,numeric,growthObjIncr-method}
\alias{growth,numeric,numeric,numeric,growthObjIncrDeclineVar-method}
\alias{growth,numeric,numeric,numeric,growthObjLogIncr-method}
\alias{growth,numeric,numeric,numeric,growthObjLogIncr.declinevar-method}
\alias{growth,numeric,numeric,numeric,growthObjTruncIncr-method}
\title{ ~~ Methods for Function \code{growth}  ~~}
\description{
 ~~ Methods for function \code{growth}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and covariates defined in data.frame. The prediction 
supplies mean sizeNext directly; the density function around this is normal. 
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov.incr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and covariates defined in data.frame. The prediction 
supplies mean incr which is added to size to obtain sizeNext; the density function around this is normal. 
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCov.logincr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and covariates defined in data.frame. The prediction 
supplies logincr which is exponentiated and added to size to obtain sizeNext; the density function around this is lognormal. 
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObj")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean sizeNext directly; the density function around this is normal. 
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjDeclineVar")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean sizeNext directly; the density function around this is normal and the sd may change with size. 
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjHossfeld")}}{
Methods to predict the probability density of \code{sizeNext} using the Hossfeld function. The prediction
supplies mean incr which is added to size to obtain sizeNext; the  density function  around this is normal. 
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjIncr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean incr which is added to size to obtain sizeNext; the density function around this is normal. 
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjIncrDeclineVar")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean incr which is added to size to obtain sizeNext; the density function around this is normal and the sd may change with size.
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjLogIncr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a  range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean logincr which is exponentiated and added to size to obtain sizeNext; the density function around this is log normal; and variance changes with size. 
%%  ~~describe this method here~~
}

\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjLogIncr.declinevar")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies logincr which is exponentiated and added to size to obtain sizeNext; the density function around this is lognormal; 
the sd may change with size. 
}


\item{\code{signature(size = "numeric", sizeNext = "numeric", cov = "numeric", growthObj = "growthObjTruncIncr")}}{
Methods to predict the probability density of \code{sizeNext} via linear prediction based around a 
range of transforms of the current size, and potentially one discrete covariate. The prediction 
supplies mean incr which is added to size to obtain sizeNext; the density function around this is truncated normal, truncated at zero. 

}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
