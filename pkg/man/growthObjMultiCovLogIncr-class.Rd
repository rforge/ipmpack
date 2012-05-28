\name{growthObjMultiCovLogIncr-class}
\Rdversion{1.1}
\docType{class}
\alias{growthObjMultiCovLogIncr-class}
%\alias{growth,numeric,numeric,data.frame,growthObjMultiCovLogIncr-method}

\title{Class \code{"growthObjMultiCovLogIncr"}}
\description{
A class object description
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("growthObjMultiCovLogIncr", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"lm"} ~~ }
     \item{\code{sd}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{growth}{\code{signature(size = "numeric", sizeNext = "numeric", cov = "data.frame", growthObj = "growthObjMultiCovLogIncr")}: ... }
	 }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("growthObjMultiCovLogIncr")
}
\keyword{classes}
